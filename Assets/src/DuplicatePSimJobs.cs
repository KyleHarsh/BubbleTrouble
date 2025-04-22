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

public class DuplicatePSimJobs : MonoBehaviour
{
    NativeArray<float3> normals;       // per‑particle normal accumulator
    NativeArray<float3> positions;
    NativeArray<float3> newPositions;
    NativeArray<float3> velocities;
    NativeArray<float3> prevPositions;
    NativeArray<float> densities;
    NativeArray<float> nearDensities;
    NativeArray<float> pressures;
    NativeArray<float> nearPressures;
    NativeArray<float3> dragRelativePositions;

    private int numParticles;

    private List<InternalParticle> particles;  // combination of surface particles and internal particles
    private NativeArray<InternalSpring> internalSprings;

    private int numSprings;
    private NativeArray<float3> springDisp;

    private NativeHashGrid neighborGrid;
    private NativeParallelMultiHashMap<int, int> neighborMap => neighborGrid.map;
    public float gridCellSize = 2.0f;
    private NativeArray<NeighborPair> neighborPairs;
    private int numPairs;

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

    private int lastSpringUpdateFrame = -1; // to limit updates to once per frame
    public float springUpdateThreshold = 0.1f; // minimum movement needed to recheck springs

    public float surfaceTensionCoeff = 0.5f;

    #region Rendering
    public Mesh particleMesh;
    public Material particleMaterial;
    public float particleRenderRadius = 0.1f;
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
    /*
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
        if (bubbleSurface != null && bubbleSurface.gameObject.activeSelf)
        {
            particles.AddRange(surfaceParticles);
        }
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

        normals = new NativeArray<float3>(numParticles, Allocator.Persistent);
        positions = new NativeArray<float3>(numParticles, Allocator.Persistent);
        newPositions = new NativeArray<float3>(numParticles, Allocator.Persistent);
        velocities = new NativeArray<float3>(numParticles, Allocator.Persistent);
        prevPositions = new NativeArray<float3>(numParticles, Allocator.Persistent);
        densities = new NativeArray<float>(numParticles, Allocator.Persistent);
        nearDensities = new NativeArray<float>(numParticles, Allocator.Persistent);
        pressures = new NativeArray<float>(numParticles, Allocator.Persistent);
        nearPressures = new NativeArray<float>(numParticles, Allocator.Persistent);
        dragRelativePositions = new NativeArray<float3>(numParticles, Allocator.Persistent);

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

        // (in SpawnParticles or right before your first substep)
        var tempList = new List<NeighborPair>();
        var iter = default(NativeParallelMultiHashMapIterator<int>);
        foreach (var kv in neighborMap) // this enumerates all (cellKey, particleIndex)
        {
            int cellKey = kv.Key;
            if (!neighborMap.TryGetFirstValue(cellKey, out int a, out iter)) continue;
            do
            {
                // for each pair (a,b) in the same cell
                if (a < kv.Value)
                    tempList.Add(new NeighborPair(a, kv.Value));
            } while (neighborMap.TryGetNextValue(out int b, ref iter));
        }
        neighborPairs = new NativeArray<NeighborPair>(tempList.Count, Allocator.Persistent);
        for (int k = 0; k < tempList.Count; k++)
            neighborPairs[k] = tempList[k];
        numPairs = neighborPairs.Length;


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
                if (j > i)
                {
                    springList.Add(new InternalSpring { particleA = i, particleB = j, restLength = r });
                }
                if (r < influenceRadius)
                {

                    float q = r / influenceRadius;
                    float w = 1f - q;
                    rho += w * w;
                }

            }
            sumRho += rho;
        }
        neighbors.Dispose();

        n0 = (sumRho / numParticles) * 2.5f;

        // upload to a native array
        internalSprings = new NativeArray<InternalSpring>(springList.Count, Allocator.Persistent);
        for (int i = 0; i < springList.Count; i++)
            internalSprings[i] = springList[i];

        numSprings = internalSprings.Length;

        springDisp = new NativeArray<float3>(numParticles, Allocator.Persistent);
    }

    private void Update()
    {

        if (particles != null)
        {
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


            { // SPH fluid simulation using substeps
                float subDeltaTime = Time.deltaTime / substeps;
                Vector3 windAcceleration = windDirection * windForce;
                Vector3 mousePos = Vector3.zero;
                Vector3 mouseWorld = Vector3.zero;

                // mouse dragging detection
                {
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
                }


                // 1) preStep: scale down velocities
                new PreStepJob
                {
                    invSS = 1f / substeps,
                    velocities = velocities
                }.Schedule(numParticles, 64).Complete();

                new ExternalForcesJob
                {
                    dt = subDeltaTime,
                    gravity = gravityForce,
                    windAccel = windAcceleration,
                    dragCoeff = dragCoefficient,
                    mousePos = mouseWorld,
                    mouseRadius = mouseDragRadius,
                    mouseStrength = mouseDragStrength,
                    positions = positions,
                    velocities = velocities,
                    prevPos = prevPositions,
                    dragRelPos = dragRelativePositions
                }.Schedule(numParticles, 64).Complete();

                for (int step = 0; step < substeps; step++)
                {
                    // a) zero‐out accumulators
                    for (int i = 0; i < normals.Length; i++) normals[i] = float3.zero;
                    for (int i = 0; i < densities.Length; i++) densities[i] = 0f;
                    for (int i = 0; i < nearDensities.Length; i++) nearDensities[i] = 0f;

                    // b) density+normal pass
                    new DensityJob
                    {
                        neigh = neighborPairs,
                        positions = positions,
                        densities = densities,
                        nearDensities = nearDensities,
                        normals = normals,
                        h = influenceRadius
                    }.Schedule(numPairs, 64).Complete();

                    // c) compute pressures on CPU
                    for (int i = 0; i < numParticles; i++)
                    {
                        pressures[i] = k_p * (densities[i] - n0);
                        nearPressures[i] = k_p2 * (nearDensities[i]);
                    }

                    // d) double‐density relaxation (PBD style)
                    new PressureRelax
                    {
                        neigh = neighborPairs,
                        positionsIn = positions,
                        positionsOut = newPositions,
                        pressures = pressures,
                        nearPressures = nearPressures,
                        dt = subDeltaTime,
                        h = influenceRadius
                    }.Schedule(numPairs, 64).Complete();

                    // swap buffers
                    var tmp = positions; positions = newPositions; newPositions = tmp;

                    // e) optional surface‐tension / adhesion
                    new SurfaceTensionPass
                    {
                        normals = normals,
                        positions = positions,
                        dt = subDeltaTime,
                        c = surfaceTensionCoeff
                    }.Schedule(numParticles, 64).Complete();

                }

                // 3) postStep: restore velocity scale
                new PostStepJob
                {
                    ss = substeps,
                    velocities = velocities
                }.Schedule(numParticles, 64).Complete();

                // 4) velocity update
                new VelocityUpdateJob
                {
                    dt = Time.deltaTime,
                    positions = positions,
                    prevPos = prevPositions,
                    velocities = velocities
                }.Schedule(numParticles, 64).Complete();

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
        normals.Dispose();
        positions.Dispose();
        newPositions.Dispose();
        velocities.Dispose();
        prevPositions.Dispose();
        densities.Dispose();
        pressures.Dispose();
        nearPressures.Dispose();
        dragRelativePositions.Dispose();
        internalSprings.Dispose();
    }*/

}

/*
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
struct PreStepJob : IJobParallelFor
{
    public float invSS;
    public NativeArray<float3> velocities;
    public void Execute(int i)
    {
        velocities[i] *= invSS;
    }
}

[BurstCompile]
struct ExternalForcesJob : IJobParallelFor
{
    public float dt;
    public float3 gravity;      // e.g. (0,-9.81,0)
    public float3 windAccel;    // your windDirection*windForce
    public float dragCoeff;    // simple air drag
    public float3 mousePos;     // world‐space mouse attractor
    public float mouseRadius;
    public float mouseStrength;

    public NativeArray<float3> positions;
    public NativeArray<float3> velocities;
    public NativeArray<float3> prevPos;
    [ReadOnly] public NativeArray<float3> dragRelPos;

    public void Execute(int i)
    {
        // save old position for velocity update later
        prevPos[i] = positions[i];

        // integrate gravity & wind
        float3 v = velocities[i];
        v += (gravity + windAccel) * dt;

        // simple air resistance
        v *= math.max(0f, 1f - dragCoeff * dt);

        // optional mouse‐drag / adhesion
        float3 toTarget = (mousePos + dragRelPos[i]) - positions[i];
        float dist = math.length(toTarget);
        if (dist < mouseRadius && dist > 1e-5f)
            v += math.normalize(toTarget) * (mouseStrength * dt);

        // write back
        velocities[i] = v;
        positions[i] += v * dt;    // predict new position
    }
}


[BurstCompile]
struct DensityJob : IJobParallelFor
{
    [ReadOnly] public NativeArray<NeighborPair> neigh;
    [ReadOnly] public NativeArray<float3> positions;
    public NativeArray<float> densities;      // ∑ w²
    public NativeArray<float> nearDensities;  // ∑ w³
    public NativeArray<float3> normals;       // ∑ w * r̂
    public float h;
    public void Execute(int k)
    {
        var np = neigh[k];
        float3 xi = positions[np.i], xj = positions[np.j];
        float3 d = xj - xi;
        float r = math.length(d);
        if (r >= h || r <= 0f) return;
        float q = r / h;
        float w = 1f - q;
        float w2 = w * w, w3 = w2 * w;
        float3 rhat = d / r;

        // **accumulate**  (if running in parallel you’ll need atomics or a two‐phase stream)
        densities[np.i] += w2;
        densities[np.j] += w2;
        nearDensities[np.i] += w3;
        nearDensities[np.j] += w3;
        normals[np.i] += w * rhat;
        normals[np.j] -= w * rhat;
    }
}


[BurstCompile]
struct PressureRelax : IJobParallelFor
{
    [ReadOnly] public NativeArray<NeighborPair> neigh;
    [ReadOnly] public NativeArray<float3> positionsIn;
    [WriteOnly] public NativeArray<float3> positionsOut;
    [ReadOnly] public NativeArray<float> pressures;
    [ReadOnly] public NativeArray<float> nearPressures;
    public float dt, h;
    public void Execute(int k)
    {
        var np = neigh[k];
        float3 xi = positionsIn[np.i], xj = positionsIn[np.j];
        float3 d = xj - xi; float r = math.length(d);
        if (r >= h || r <= 0f) return;
        float q = r / h;
        float w0 = 1f - q, w1 = w0 * w0;
        float Pi = pressures[np.i], Pj = pressures[np.j];
        float Nip = nearPressures[np.i], Npj = nearPressures[np.j];
        float coeff = dt * dt * ((Pi + Pj) * w0 + (Nip + Npj) * w1);
        float3 dir = d / r;
        float3 D = coeff * dir;
        // split
        positionsOut[np.i] -= 0.5f * D;
        positionsOut[np.j] += 0.5f * D;
    }
}

[BurstCompile]
struct SurfaceTensionPass : IJobParallelFor
{
    [ReadOnly] public NativeArray<float3> normals;
    public NativeArray<float3> positions;
    public float dt, c;  // surface tension coeff
    public void Execute(int i)
    {
        float3 n = normals[i];
        float len = math.length(n);
        if (len < 1e-6f) return;
        float curvature = -len; // approximate mean curvature
        float3 D = dt * dt * c * curvature * (n / len);
        positions[i] += D;
    }
}


[BurstCompile]
struct VelocityUpdateJob : IJobParallelFor
{
    public float dt;
    [ReadOnly] public NativeArray<float3> positions, prevPos;
    public NativeArray<float3> velocities;

    public void Execute(int i)
    {
        velocities[i] = (positions[i] - prevPos[i]) / dt;
    }
}

[BurstCompile]
struct PostStepJob : IJobParallelFor
{
    public float ss;
    public NativeArray<float3> velocities;
    public void Execute(int i)
    {
        velocities[i] *= ss;
    }
}


[BurstCompile]
struct SpringJob : IJob
{
    public float dt;
    public float k_spring;
    public float h;                    // influence radius
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

            // Hooke force + damper
            float Fs = k_spring * (s.restLength - r);

            float3 dp = dt * dt * (1f - s.restLength / h) * Fs * dir;

            positionsOut[s.particleA] -= dp * 0.5f;
            positionsOut[s.particleB] += dp * 0.5f;
        }

    }
}

[BurstCompile]
struct AccumulateJob : IJobParallelFor
{
    public float cellSize;
    [ReadOnly] public NativeArray<float3> positions;
    public NativeArray<float3> normals;
    public NativeArray<float> densities;
    [ReadOnly] public NativeParallelMultiHashMap<int, int> map;

    public void Execute(int i)
    {
        float3 pi = positions[i];
        int cx = (int)math.floor(pi.x / cellSize),
            cy = (int)math.floor(pi.y / cellSize),
            cz = (int)math.floor(pi.z / cellSize);

        float3 Ni = float3.zero;
        float Di = 0f;

        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                {
                    int key = HashCell(cx + dx, cy + dy, cz + dz);
                    if (map.TryGetFirstValue(key, out int j, out var it))
                    {
                        do
                        {
                            if (j == i) continue;
                            float3 pj = positions[j];
                            float3 rij = pi - pj;
                            float r2 = math.dot(rij, rij);
                            if (r2 >= cellSize * cellSize) continue;
                            float r = math.sqrt(r2);
                            float w = math.max(0f, 1f - r / cellSize);
                            float3 n = rij / (r + 1e-6f);        // avoid div0
                            Ni += w * n;                         // accumulate normal
                            Di += w * w;                           // accumulate density
                        } while (map.TryGetNextValue(out j, ref it));
                    }
                }

        normals[i] = Ni;
        densities[i] = Di;
    }

    static int HashCell(int x, int y, int z)
    {
        unchecked
        {
            int h = x * 73856093; h ^= y * 19349663; h ^= z * 83492791;
            return h;
        }
    }
}

[BurstCompile]
struct ComputePressureJob : IJobParallelFor
{
    public float n0, k_p, k_p2;
    public float cellSize;
    [ReadOnly] public NativeArray<float3> normals;
    [ReadOnly] public NativeArray<float> densities;
    public NativeArray<float> pressures;
    public NativeArray<float> nearPressures;

    public void Execute(int i)
    {
        float D = densities[i];
        pressures[i] = k_p * (D - n0);
        nearPressures[i] = k_p2 * math.pow(math.length(normals[i]) / cellSize, 3);
        // note: Drops uses w³ kernel for nearPressure; this approximates it
    }
}

[BurstCompile]
struct PressureImpulseJob : IJobParallelFor
{
    public float dt, cellSize;
    [ReadOnly] public NativeArray<float3> positions;
    public NativeArray<float3> velocities;
    [ReadOnly] public NativeArray<float> pressures;
    [ReadOnly] public NativeArray<float> nearPressures;
    [ReadOnly] public NativeParallelMultiHashMap<int, int> map;

    public void Execute(int i)
    {
        float3 pi = positions[i];
        float3 vi = velocities[i];

        int cx = (int)math.floor(pi.x / cellSize),
            cy = (int)math.floor(pi.y / cellSize),
            cz = (int)math.floor(pi.z / cellSize);

        float3 Δv = float3.zero;

        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                {
                    int key = HashCell(cx + dx, cy + dy, cz + dz);
                    if (map.TryGetFirstValue(key, out int j, out var it))
                    {
                        do
                        {
                            if (j == i) continue;
                            float3 pj = positions[j];
                            float3 rij = pi - pj;
                            float r2 = math.dot(rij, rij);
                            if (r2 >= cellSize * cellSize) continue;
                            float r = math.sqrt(r2);
                            float w = math.max(0f, 1f - r / cellSize);
                            float Pij = (pressures[i] + pressures[j] + nearPressures[i] + nearPressures[j]);
                            float3 n = rij / (r + 1e-6f);
                            Δv += 0.5f * dt * dt * w * Pij * n;
                        } while (map.TryGetNextValue(out j, ref it));
                    }
                }

        velocities[i] = vi - Δv;  // subtract because we did pi-pj
    }

    static int HashCell(int x, int y, int z)
    {
        unchecked
        {
            int h = x * 73856093; h ^= y * 19349663; h ^= z * 83492791;
            return h;
        }
    }
}


#endregion

*/