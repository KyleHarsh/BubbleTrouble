using NUnit.Framework;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.VisualScripting;
using UnityEngine;
using static UnityEditor.PlayerSettings;
using static UnityEngine.ParticleSystem;

public class ParticleSimulatorJobs : MonoBehaviour
{
    NativeArray<float3> positions;  // current position
    NativeArray<float3> newPositions;  // buffer position for multi threading
    NativeArray<float3> velocities;    // used for initial external force prediction
    NativeArray<float3> prevPositions;  // position values at start of substep
    NativeArray<float> densities;
    NativeArray<float> pressures;
    NativeArray<float> nearPressures;
    NativeArray<float3> dragRelativePositions;
    NativeArray<float3> normals;       // ∇C
    NativeArray<float> colorField;    // C

    private int numParticles;

    private List<InternalParticle> particles;  // combination of surface particles and internal particles
    private NativeArray<InternalSpring> internalSprings;

    private int numSprings;
    private NativeArray<float3> springDisp;

    private NativeHashGrid neighborGrid;
    private NativeParallelMultiHashMap<int, int> neighborMap => neighborGrid.map;
    public float gridCellSize = 2.0f;
    public BubbleParticleSpawner spawner;

    public Transform cameraTransform;
    public Transform mouseScreenTransform;

    private bool spawned = false;

    public int substeps = 4;

    public float k_p = 0.004f;   // Pressure stiffness coefficient
    public float k_p2 = 0.01f;   // Near-pressure coefficient
    public float n0 = 10f;       // Rest (target) density
    public float influenceRadius;
    public int pressureIterations = 3;

    public float gravityForce = 9.8f;
    public float particleMass = 1.0f;

    public float springStiffness = 0.3f;  // k_spring, controls how strong the spring force is.
    public float plasticity = 0.2f;
    public float yieldRatio = 0.1f;
    public float dampingCoeff = 0.2f;

    private int lastSpringUpdateFrame = -1; // to limit updates to once per frame
    public float springUpdateThreshold = 0.1f; // minimum movement needed to recheck springs

    public float surfaceTensionCoeff = 0.5f;

    #region Rendering
    public Mesh particleMesh;
    public Material particleMaterial;
    public float particleRenderRadius = 0.1f;
    public bool renderParticles = true;
    #endregion

    private List<Matrix4x4> transformMatrices = new List<Matrix4x4>();

    private float distanceToCamera;
    // Define a radius for influence for dragging
    public float mouseDragRadius = 1.0f;
    public float mouseDragStrength = 1.0f;
    [SerializeField] private bool dragging = false;
    public float dragMomentumDuration = 1f;
    [SerializeField] private float dragTimer = 0.0f;

    public float dragCoefficient = 0.1f; // Tweak this parameter

    public Vector3 windDirection;
    public float windForce;

    [SerializeField] private Bubble_m bubbleSurface;
    private List<SurfaceParticle> surfaceParticles;

    // Start is called before the first frame update
    void Start()
    {
        if (bubbleSurface != null && bubbleSurface.gameObject.activeSelf)
        {
            bubbleSurface.transform.position = spawner.transform.position;
            surfaceParticles = bubbleSurface.surfaceParticles;
        }
            
        distanceToCamera = Vector3.Distance(cameraTransform.position, mouseScreenTransform.position);
        SpawnParticles();
    }

    public void SpawnParticles()
    {
        particles = spawner.SpawnInternalParticles();
        if (bubbleSurface != null && bubbleSurface.gameObject.activeSelf) particles.AddRange(surfaceParticles);
        numParticles = particles.Count;

        float invDt2 = substeps * substeps; // = 16
        k_p *= invDt2;   // e.g. from 0.004 → 0.064
        k_p2 *= invDt2;   // e.g. from 0.01  → 0.16

        // 2) Compute Poisson spacing & smoothing length (h = 2×spacing)
        float sphereVolume = (4f / 3f) * Mathf.PI * Mathf.Pow(spawner.bubbleRadius, 3);
        float avgVolPerPt = sphereVolume / numParticles;
        float spacing = Mathf.Pow(avgVolPerPt, 1f / 3f);
        influenceRadius = spacing * 2f;

        // Use same h for the grid
        neighborGrid = new NativeHashGrid(numParticles * 2, influenceRadius, Allocator.Persistent);

        positions = new NativeArray<float3>(numParticles, Allocator.Persistent);
        newPositions = new NativeArray<float3>(numParticles, Allocator.Persistent);
        velocities = new NativeArray<float3>(numParticles, Allocator.Persistent);
        prevPositions = new NativeArray<float3>(numParticles, Allocator.Persistent);
        densities = new NativeArray<float>(numParticles, Allocator.Persistent);
        pressures = new NativeArray<float>(numParticles, Allocator.Persistent);
        nearPressures = new NativeArray<float>(numParticles, Allocator.Persistent);
        dragRelativePositions = new NativeArray<float3>(numParticles, Allocator.Persistent);
        normals = new NativeArray<float3>(numParticles, Allocator.Persistent);
        colorField = new NativeArray<float>(numParticles, Allocator.Persistent);

        // 4) Copy initial data into arrays
        for (int i = 0; i < numParticles; i++)
        {
            var p = particles[i];
            positions[i] = p.position;
            velocities[i] = p.velocity;
            prevPositions[i] = p.position;
        }

        // 5) Build the grid once, sample rest density n0
        neighborGrid.Clear();
        for (int i = 0; i < numParticles; i++)
            neighborGrid.Insert(positions[i], i);

        var springList = new List<InternalSpring>();
        var neighbors = new NativeList<int>(Allocator.Temp);
        float sumRho = 0f;
        for (int i = 0; i < numParticles; i++)
        {
            neighbors.Clear();
            float rho = 1f;
            neighborGrid.GetNeighbors(positions[i], ref neighbors);

            foreach (int j in neighbors)
            {
                float r = math.distance(positions[i], positions[j]);

                if (r < influenceRadius)
                {
                    springList.Add(new InternalSpring { particleA = i, particleB = j, restLength = r });

                    float q = r / influenceRadius;
                    float w = 1f - q;
                    rho += w * w;
                }

            }
            sumRho += rho;
        }
        neighbors.Dispose();

        n0 = (sumRho / numParticles)*2.5f;

        // upload to a native array
        internalSprings = new NativeArray<InternalSpring>(springList.Count, Allocator.Persistent);
        for (int i = 0; i < springList.Count; i++)
            internalSprings[i] = springList[i];

        numSprings = internalSprings.Length;

        //springDisp = new NativeArray<float3>(numParticles, Allocator.Persistent);
    }

    private void Update()
    {

        if (particles != null)
        {
            if(renderParticles)
            { // rendering

                // Prepare instance matrices for each particle.
                int count = particles.Count;
                Matrix4x4[] matrices = new Matrix4x4[count];
                for (int i = 0; i < count; i++)
                {
                    // transformation matrix for each particle
                    matrices[i] = Matrix4x4.TRS(particles[i].position, Quaternion.identity, Vector3.one * particleRenderRadius);
                }

                // Render all particles in one call
                Graphics.DrawMeshInstanced(particleMesh, 0, particleMaterial, matrices);
            }
            //Debug.Log($"n₀ = {n0:F2}, ρ[0] = {densities[0]:F2},  ρₘᵢₙ = {densities.Min():F2}, ρₘₐₓ = {densities.Max():F2}");


            { // SPH fluid simulation using substeps
                float subDeltaTime = Time.deltaTime / substeps;
                Vector3 windAcceleration = windDirection * windForce;
                Vector3 mousePos = Vector3.zero;
                Vector3 mouseWorld = Vector3.zero;

                if (!dragging && dragTimer > 0)
                {
                    dragTimer -= Time.deltaTime;
                }

                // convert screen→world at some plane or depth
                if (Input.GetMouseButtonDown(0) && !dragging)
                {
                    dragging = true;
                    for (int i = 0; i < particles.Count; i++)
                    {
                        Vector3 clickPos = Input.mousePosition;
                        clickPos.z = cameraTransform.position.z * -1f;
                        clickPos = Camera.main.ScreenToWorldPoint(clickPos);
                        if (Vector3.Distance(positions[i], clickPos) < mouseDragRadius)
                        {
                            Vector3 vecPos = new Vector3(positions[i].x, positions[i].y, positions[i].z);
                            dragRelativePositions[i] = vecPos - clickPos;
                        }
                        else
                        {
                            dragRelativePositions[i] = float3.zero;
                        }
                    }
                }
                if (Input.GetMouseButton(0) && dragging)
                {
                    mousePos = Input.mousePosition;
                    mousePos.z = cameraTransform.position.z * -1f; // or set to your fluid depth
                    mouseWorld = Camera.main.ScreenToWorldPoint(mousePos);
                }
                if (Input.GetMouseButtonUp(0) && dragging)
                {
                    dragging = false;
                    dragTimer = dragMomentumDuration;
                }

                NativeArray<float3> springCorrected = new NativeArray<float3>(numParticles, Allocator.TempJob);

                // B) Rebuild grid with new positions
                neighborGrid.Clear();
                for (int i = 0; i < numParticles; i++)
                    neighborGrid.Insert(positions[i], i);

                for (int step = 0; step < substeps; step++)
                {
                    // A) Predict
                    var pj = new PredictJob
                    {
                        gravity = new float3(0, -gravityForce, 0),
                        dt = subDeltaTime,
                        dragCoeff = dragCoefficient,
                        windAccel = windAcceleration,
                        mousePos = mouseWorld,
                        mouseRadius = mouseDragRadius,
                        mouseStrength = (dragging) ? mouseDragStrength : ((dragTimer > 0)? mouseDragStrength*(dragTimer/dragMomentumDuration) : 0f),
                        positions = positions,
                        velocities = velocities,
                        prevPos = prevPositions,
                        dragRelativePositions = dragRelativePositions
                    };
                    JobHandle h1 = pj.Schedule(numParticles, 64);
                    h1.Complete(); // we need positions[] updated before building the grid

                    neighborGrid.Clear();
                    for (int i = 0; i < numParticles; i++)
                        neighborGrid.Insert(positions[i], i);

                    // rebuild spring list
                    List<InternalSpring> springList = new List<InternalSpring>();
                    var nbrs = new NativeList<int>(Allocator.Temp);
                    for (int i = 0; i < numParticles; i++)
                    {
                        nbrs.Clear();
                        neighborGrid.GetNeighbors(positions[i], ref nbrs);
                        foreach (int j in nbrs)
                        {
                            float r = math.distance(positions[i], positions[j]);
                            if(r < influenceRadius)
                            springList.Add(new InternalSpring
                            {
                                particleA = i,
                                particleB = j,
                                restLength = r
                            });
                        }
        
                    }
                    nbrs.Dispose();

                    // plasticity & prune
                    for (int s = 0; s < springList.Count; s++)
                    {
                        var sp = springList[s];
                        float r = math.distance(positions[sp.particleA], positions[sp.particleB]);
                        float tol = yieldRatio * sp.restLength;
                        if (r > sp.restLength + tol)
                            sp.restLength += plasticity * (r - sp.restLength - tol) * subDeltaTime;
                        else if (r < sp.restLength - tol)
                            sp.restLength -= plasticity * (sp.restLength - tol - r) * subDeltaTime;

                        if (sp.restLength > influenceRadius)
                        {
                            springList.RemoveAt(s--);
                            continue;
                        }

                        springList[s] = sp;
                    }

                    // copy back into your NativeArray
                    internalSprings.Dispose();
                    internalSprings = new NativeArray<InternalSpring>(springList.Count, Allocator.TempJob);
                    for (int i = 0; i < springList.Count; i++)
                        internalSprings[i] = springList[i];

                    // E) Spring corrections 
                    // 1. copy current positions into the scratch buffer

                    springCorrected.CopyFrom(positions);

                    var springJob = new SpringJob()
                    {
                        dt = subDeltaTime,
                        h = influenceRadius,
                        k_spring = springStiffness,
                        dampForce = dampingCoeff,
                        springs = internalSprings,
                        positionsIn = positions,
                        positionsOut = springCorrected
                    };
                    JobHandle h4 = springJob.Schedule(h1);
                    h4.Complete();

                    positions.CopyFrom(springCorrected);

                    neighborGrid.Clear();
                    for (int i = 0; i < numParticles; i++)
                        neighborGrid.Insert(positions[i], i);

                    // … after spring corrections …
                    springCorrected.CopyFrom(positions);

                    var stJob = new SurfaceTensionJob
                    {
                        h = influenceRadius,
                        sigma = surfaceTensionCoeff,
                        dt = subDeltaTime,
                        positions = positions,
                        map = neighborMap,
                        colorField = colorField,
                        normals = normals,
                        outPos = springCorrected
                    };
                    JobHandle stj = stJob.Schedule(numParticles, 64, h4);
                    stj.Complete();

                    // copy back
                    positions.CopyFrom(springCorrected);

                    neighborGrid.Clear();
                    for (int i = 0; i < numParticles; i++)
                        neighborGrid.Insert(positions[i], i);

                    // C) Density
                    var dj = new DensityJob
                    {
                        cellSize = influenceRadius,
                        positions = positions,
                        densities = densities,
                        pressures = pressures,
                        nearPressures = nearPressures,
                        n0 = n0,
                        k_p = k_p,
                        k_p2 = k_p2,
                        map = neighborMap
                    };
                    JobHandle h2 = dj.Schedule(numParticles, 64, stj);
                    h2.Complete();

                    // 4) *pressure solve iterations*
                    JobHandle h3;
                    for (int iter = 0; iter < pressureIterations; iter++)
                    {
                        var disp = new DisplacementJob
                        {
                            cellSize = influenceRadius,
                            dt = subDeltaTime,
                            positionsIn = positions,
                            positionsOut = newPositions,
                            pressures = pressures,
                            nearPressures = nearPressures,
                            map = neighborMap
                        };
                        h3 = disp.Schedule(numParticles, 64, h2);
                        h3.Complete();

                        // swap buffers so next iteration sees the corrected positions
                        var temp = positions; positions = newPositions; newPositions = temp;

                        // recompute densities & pressures on the moved positions
                        dj = new DensityJob
                        {
                            cellSize = influenceRadius,
                            positions = positions,
                            densities = densities,
                            pressures = pressures,
                            nearPressures = nearPressures,
                            n0 = n0,
                            k_p = k_p,
                            k_p2 = k_p2,
                            map = neighborMap
                        };
                        h2 = dj.Schedule(numParticles, 64);
                        h2.Complete();
                    }
                    var vu = new VelocityUpdateJob
                    {
                        dt = subDeltaTime,
                        positions = newPositions,
                        prevPositions = prevPositions,
                        velocities = velocities
                    };
                    JobHandle h5 = vu.Schedule(numParticles, 64, h2);
                    h5.Complete();   
                }

                springCorrected.Dispose();
                
                // 2) Copy back for rendering
                for (int i = 0; i < numParticles; i++)
                {
                    particles[i].position = positions[i];
                    particles[i].velocity = velocities[i];
                }


            }
        }

    }

    #region External Forces

    void ApplyMouseDragForce(float deltaTime)
    {
        if (Input.GetMouseButton(0))
        {  // If the mouse button is pressed
            Debug.Log("dragging");
            // Convert the mouse position into a world coordinate.
            Vector3 mouseScreenPos = Input.mousePosition;
            mouseScreenPos.z = distanceToCamera;  // Set appropriate depth.
            Vector3 mouseWorldPos = Camera.main.ScreenToWorldPoint(mouseScreenPos);

            // Iterate over all particles.
            foreach (var particle in particles)
            {
                float distance = Vector3.Distance(particle.position, mouseWorldPos);
                if (distance < mouseDragRadius)
                {
                    // Compute a simple drag force.
                    Vector3 force = mouseDragStrength * (mouseWorldPos - particle.position);
                    // You can either add this force to the particle’s velocity or adjust its position directly.
                    particle.position += force * deltaTime * deltaTime;
                }

            }
        }
    }
    #endregion

    private void OnDrawGizmos()
    {
        if (Input.GetMouseButton(0))
        {
            // Convert the mouse position into a world coordinate.
            Vector3 mouseScreenPos = Input.mousePosition;
            mouseScreenPos.z = distanceToCamera;  // Set appropriate depth.
            Vector3 mouseWorldPos = Camera.main.ScreenToWorldPoint(mouseScreenPos);

            Gizmos.color = Color.cyan.WithAlpha(0.5f);
            Gizmos.DrawSphere(mouseWorldPos, mouseDragRadius);
        }
    }

    private void OnDestroy()
    {
        positions.Dispose();
        newPositions.Dispose();
        velocities.Dispose();
        prevPositions.Dispose();
        densities.Dispose();
        pressures.Dispose();
        nearPressures.Dispose();
        dragRelativePositions.Dispose();
        internalSprings.Dispose();
        springDisp.Dispose();
        normals.Dispose();
        colorField.Dispose();
    }

}


#region Job Structs
[BurstCompile]
struct PredictJob : IJobParallelFor
{
    public float3 gravity;
    public float dt;
    public float3 windAccel;
    public float dragCoeff;
    public float3 mousePos;
    public float mouseRadius;
    public float mouseStrength;

    public NativeArray<float3> positions;
    public NativeArray<float3> velocities;
    public NativeArray<float3> prevPos;
    public NativeArray<float3> dragRelativePositions;

    public void Execute(int i)
    {
        // 1) save old
        prevPos[i] = positions[i];

        float3 v = velocities[i];

        v += (gravity + windAccel) * dt;

        v *= (1f - dragCoeff) * dt;

        // new changes, particle drag target is set point away from mouse position when click start (within mouseRadius)

        float3 dragTarget = mousePos + dragRelativePositions[i];
        float3 diff = dragTarget - positions[i];
        float d = math.length(diff);
        if (d < mouseRadius && d > 1e-5f)
        {
            float3 dir = diff / d;
           
            v += dir * (mouseStrength * dt);
        }

        velocities[i] = v;
        // positions array is your “current” buffer so the next jobs see it
        positions[i] += v * dt;
    }
}

[BurstCompile]
struct DensityJob : IJobParallelFor
{
    public float cellSize;
    [ReadOnly] public NativeArray<float3> positions;
    [WriteOnly] public NativeArray<float> densities, pressures, nearPressures;
    public float n0, k_p, k_p2;
    [ReadOnly] public NativeParallelMultiHashMap<int, int> map;

    public void Execute(int i)
    {
        float3 pi = positions[i];
        int cx = (int)math.floor(pi.x / cellSize);
        int cy = (int)math.floor(pi.y / cellSize);
        int cz = (int)math.floor(pi.z / cellSize);

        // in DensityJob.Execute:
        float sum = 1f;       // ← start at 1 for the self‐density
        float sumNear = 1f;   // ← ditto
        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                {
                    int key = HashCell(cx + dx, cy + dy, cz + dz);
                    if (map.TryGetFirstValue(key, out int j, out var iter))
                    {
                        do
                        {
                            float3 pj = positions[j];
                            float r = math.distance(pi, pj);
                            if (r < cellSize)
                            {
                                float q = r / cellSize;
                                float w = 1f - q;
                                sum += w * w;
                                sumNear += w * w * w;
                            }
                        } while (map.TryGetNextValue(out j, ref iter));
                    }
                }

        densities[i] = sum;
        pressures[i] = k_p * (sum - n0);
        nearPressures[i] = k_p2 * sumNear;
    }

    static int HashCell(int x, int y, int z)
    {
        unchecked
        {
            int h = x * 73856093;
            h ^= y * 19349663;
            h ^= z * 83492791;
            return h;
        }
    }
}

[BurstCompile]
struct DisplacementJob : IJobParallelFor
{
    public float cellSize, dt;
    [ReadOnly] public NativeArray<float3> positionsIn;
    [WriteOnly] public NativeArray<float3> positionsOut;
    [ReadOnly] public NativeArray<float> pressures, nearPressures;
    [ReadOnly] public NativeParallelMultiHashMap<int, int> map;

    public void Execute(int i)
    {
        float3 pi = positionsIn[i];
        int cx = (int)math.floor(pi.x / cellSize);
        int cy = (int)math.floor(pi.y / cellSize);
        int cz = (int)math.floor(pi.z / cellSize);

        float3 disp = float3.zero;
        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                {
                    int key = HashCell(cx + dx, cy + dy, cz + dz);
                    if (map.TryGetFirstValue(key, out int j, out var iter))
                    {
                        do
                        {
                            if (j == i) continue;
                            float3 pj = positionsIn[j];
                            float r = math.distance(pi, pj);
                            if (r <= 0f || r >= cellSize) continue;
                            float q = r / cellSize;
                            float w0 = 1f - q;
                            float w1 = w0 * w0;
                            float P = pressures[i] + pressures[j];
                            float NP = nearPressures[i] + nearPressures[j];
                            float coeff = -(dt*dt) * (P * w0 + NP * w1);
                            float3 dir = (pj - pi) / r;
                            disp += coeff * dir;
                        } while (map.TryGetNextValue(out j, ref iter));
                    }
                }
        positionsOut[i] = pi + disp * 0.5f;
    }

    static int HashCell(int x, int y, int z) { unchecked { int h = x * 73856093; h ^= y * 19349663; h ^= z * 83492791; return h; } }
}

[BurstCompile]
struct VelocityUpdateJob : IJobParallelFor
{
    public float dt;
    [ReadOnly] public NativeArray<float3> positions, prevPositions;
    public NativeArray<float3> velocities;

    public void Execute(int i)
    {
        velocities[i] = (positions[i] - prevPositions[i]) / dt;
    }
}

[BurstCompile]
struct SpringJob : IJob
{
    public float dt;
    public float k_spring;
    public float h;                    // influence radius
    public float dampForce;
    public NativeArray<InternalSpring> springs;
    public NativeArray<float3> positionsIn;
    public NativeArray<float3> positionsOut;

    public void Execute()
    {
        for (int si = 0; si < springs.Length; si++)
        {
            var s = springs[si];
            float3 pa = positionsIn[s.particleA];
            float3 pb = positionsIn[s.particleB];

            float3 dir = pb - pa;
            float r = math.length(dir);
            if (r <= 1e-6f || r > h) continue;
            dir /= r;

            // Hooke force 
            float Fs = k_spring * (s.restLength - r);
            float Fd = k_spring * dampForce;


            float3 dp = dt * dt * (1f - s.restLength / h) * ((Fs * dir) - (Fd * -1f * dir));

            positionsOut[s.particleA] -= dp * 0.5f;
            positionsOut[s.particleB] += dp * 0.5f;
        }
        
    }
}

[BurstCompile]
struct SurfaceTensionJob : IJobParallelFor
{
    public float h;          // influence radius
    public float sigma;      // surface‑tension coefficient
    public float dt;         // timestep
    [ReadOnly] public NativeArray<float3> positions;
    [ReadOnly] public NativeParallelMultiHashMap<int, int> map;

    public NativeArray<float> colorField;   // C_i
    public NativeArray<float3> normals;      // ∇C_i
    public NativeArray<float3> outPos;       // scratch positionsOut

    static int HashCell(int x, int y, int z)
    {
        unchecked
        {
            int h1 = x * 73856093; h1 ^= y * 19349663; h1 ^= z * 83492791; return h1;
        }
    }

    public void Execute(int i)
    {
        float3 pi = positions[i];
        int cx = (int)math.floor(pi.x / h),
            cy = (int)math.floor(pi.y / h),
            cz = (int)math.floor(pi.z / h);

        // SPH sums for color field C and its gradient ∇C
        float Ci = 0f;
        float3 grad = float3.zero;

        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                {
                    int key = HashCell(cx + dx, cy + dy, cz + dz);
                    if (!map.TryGetFirstValue(key, out int j, out var it)) continue;
                    do
                    {
                        if (j == i) continue;
                        float3 pj = positions[j];
                        float r = math.distance(pi, pj);
                        if (r >= h) continue;
                        float q = r / h;

                        // use smoothing kernel W(q) = (1–q)^3
                        float w = (1 - q) * (1 - q) * (1 - q);
                        Ci += w;
                        // gradient of W : ∇W = –3*(1–q)^2 * (1/h) * (pj–pi)/r
                        float dw = -3f * (1 - q) * (1 - q) / (h * r);
                        grad += dw * (pj-pi);
                    } while(map.TryGetNextValue(out j, ref it));
                }

        colorField[i] = Ci + 1f;      // +1 to include “self”
        normals[i] = grad;         // ∇C

        // if we’re “pulling” positions directly (PBD), we can
        // turn F = –σ * κ * (∇C/|∇C|) into a position correction:
        float Nlen = math.length(grad);
        if (Nlen > 1e-5f)
        {
            // approximate curvature κ ≈ – ∇²C / |∇C| :
            // we can re‑use Ci – 1 as rough ∇²C since ∇²C≈∑W
            float lap = Ci;
            float3 force = sigma * lap * (grad / Nlen);
            // impulses ~ F * dt²
            outPos[i] += force * dt * dt;
        }
    }
}


#endregion
