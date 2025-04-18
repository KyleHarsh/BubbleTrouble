using NUnit.Framework;
using System.Collections;
using System.Collections.Generic;
using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.VisualScripting;
using UnityEngine;
using static UnityEngine.ParticleSystem;

public class ParticleSimulatorJobs : MonoBehaviour
{
    NativeArray<float3> positions;
    NativeArray<float3> newPositions;
    NativeArray<float3> velocities;
    NativeArray<float3> prevPositions;
    NativeArray<float> densities;
    NativeArray<float> pressures;
    NativeArray<float> nearPressures;

    private int numParticles;

    private List<InternalParticle> internalParticles;
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

    public float k_p = 0.004f;   // Pressure stiffness coefficient
    public float k_p2 = 0.01f;   // Near-pressure coefficient
    public float n0 = 10f;       // Rest (target) density
    public float influenceRadius;

    public float springStiffness = 0.3f;  // k_spring, controls how strong the spring force is.
    public float plasticity = 0.3f;       // α, the plasticity constant for rest length adjustment.
    public float yieldRatio = 0.1f;       // γ, the yield threshold ratio.

    private int lastSpringUpdateFrame = -1; // to limit updates to once per frame
    public float springUpdateThreshold = 0.1f; // minimum movement needed to recheck springs

    #region Rendering
    public Mesh particleMesh;
    public Material particleMaterial;
    public float particleRenderRadius = 0.1f;
    #endregion

    private List<Matrix4x4> transformMatrices = new List<Matrix4x4>();

    public float gravityForce = 9.8f;
    public float particleMass;

    private float distanceToCamera;
    // Define a radius for influence.
    public float mouseDragRadius = 1.0f;
    public float mouseDragStrength = 1.0f;

    public float dragCoefficient = 0.1f; // Tweak this parameter

    // Start is called before the first frame update
    void Start()
    {
        SpawnParticles();
    }

    public void SpawnParticles()
    {
        internalParticles = spawner.SpawnInternalParticles();
        numParticles = internalParticles.Count;

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

        // 4) Copy initial data into arrays
        for (int i = 0; i < numParticles; i++)
        {
            var p = internalParticles[i];
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
            float rho = 0f;
            neighborGrid.GetNeighbors(positions[i], ref neighbors);

            foreach (int j in neighbors)
            {
                float r = math.distance(positions[i], positions[j]);
                if (j > i)
                {
                    springList.Add(new InternalSpring { particleA = i, particleB = j, restLength = r });
                }
                else if (j == i) continue;
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

        n0 = (sumRho / numParticles) * 2f;

        // upload to a native array
        internalSprings = new NativeArray<InternalSpring>(springList.Count, Allocator.Persistent);
        for (int i = 0; i < springList.Count; i++)
            internalSprings[i] = springList[i];

        numSprings = internalSprings.Length;

        springDisp = new NativeArray<float3>(numParticles, Allocator.Persistent);
    }

    private void Update()
    {

        if (internalParticles != null)
        {
            { // rendering
                // Update your particle simulation here (or elsewhere) and update the particles list

                // Prepare instance matrices for each particle.
                int count = internalParticles.Count;
                Matrix4x4[] matrices = new Matrix4x4[count];
                for (int i = 0; i < count; i++)
                {
                    // transformation matrix for each particle
                    matrices[i] = Matrix4x4.TRS(internalParticles[i].position, Quaternion.identity, Vector3.one * particleRenderRadius);
                }

                // Render all particles in one call.
                // Note: There is a limit to the number of instances per call (e.g., 1023), so you might need to batch them.
                Graphics.DrawMeshInstanced(particleMesh, 0, particleMaterial, matrices);
            }

            { // SPH fluid simulation using substeps
                int substeps = 4;
                float subDeltaTime = Time.deltaTime / substeps;

                for (int step = 0; step < substeps; step++)
                {
                    // A) Predict
                    var pj = new PredictJob
                    {
                        gravity = new float3(0, -gravityForce, 0),
                        dt = subDeltaTime,
                        positions = positions,
                        velocities = velocities,
                        prevPos = prevPositions
                    };
                    JobHandle h1 = pj.Schedule(numParticles, 64);
                    h1.Complete(); // we need positions[] updated before building the grid

                    // B) Rebuild grid with new positions
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
                    JobHandle h2 = dj.Schedule(numParticles, 64, h1);

                    // D) Displacement
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
                    JobHandle h3 = disp.Schedule(numParticles, 64, h2);
                    h3.Complete();

                    for (int i = 0; i < pressures.Length; i++)
                    {
                        pressures[i] = nearPressures[i] = 0f;
                    }

                    // E) Spring corrections on the CPU
                    // ZERO OUT springDisp
                    for(int i = 0; i < springDisp.Length; i++)
                    {
                        springDisp[i] = float3.zero;
                    }

                    // after your DisplacementJob completes and writes newPositions…
                    var springJob = new SpringJob()
                    {
                        dt = subDeltaTime,
                        k_spring = springStiffness,
                        plasticity = plasticity,
                        yieldRatio = yieldRatio,
                        h = influenceRadius,
                        springs = internalSprings,
                        newPositions = newPositions
                    };
                    springJob.Schedule().Complete();

                    // F) Velocity update
                    var vu = new VelocityUpdateJob
                    {
                        dt = subDeltaTime,
                        positions = newPositions,
                        prevPos = prevPositions,
                        velocities = velocities
                    };
                    JobHandle h4 = vu.Schedule(numParticles, 64);
                    h4.Complete();
                    
                    /*
                    // small damping
                    for (int i = 0; i < numParticles; i++)
                        velocities[i] *= 0.99f;*/
                    
                    // G) Swap buffers for next substep
                    var tmp = positions; positions = newPositions; newPositions = tmp;
                }

                // 2) Copy back for rendering
                for (int i = 0; i < numParticles; i++)
                {
                    internalParticles[i].position = positions[i];
                    internalParticles[i].velocity = velocities[i];
                }

            }
        }

    }

    #region External Forces

    private void ApplyAirResistance(float dt)
    {
        foreach (InternalParticle p in internalParticles)
        {
            // Reduces the velocity proportionally to the current value.
            p.velocity *= (1.0f - dragCoefficient * dt);
        }
    }

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
            foreach (var particle in internalParticles)
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
        if (springDisp.IsCreated) springDisp.Dispose();
    }

}


#region Job Structs
[BurstCompile]
struct PredictJob : IJobParallelFor
{
    public float3 gravity;
    public float dt;
    public NativeArray<float3> positions, velocities, prevPos;

    public void Execute(int i)
    {
        prevPos[i] = positions[i];
        velocities[i] += gravity * dt;
        positions[i] += velocities[i] * dt;
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

        float sum = 0f, sumNear = 0f;
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
                            float coeff = (dt * dt) * (P * w0 + NP * w1);
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
    [ReadOnly] public NativeArray<float3> positions, prevPos;
    public NativeArray<float3> velocities;

    public void Execute(int i)
    {
        velocities[i] = (positions[i] - prevPos[i]) / dt;
    }
}

[BurstCompile]
struct SpringJob : IJob
{
    public float dt;
    public float k_spring;
    public float plasticity;
    public float yieldRatio;
    public float h;                    // influence radius
    public NativeArray<InternalSpring> springs;  // your NativeArray<Spring>
    public NativeArray<float3> newPositions;

    public void Execute()
    {
        for (int si = 0; si < springs.Length; si++)
        {
            var s = springs[si];

            float3 pa = newPositions[s.particleA];
            float3 pb = newPositions[s.particleB];
            float r = math.distance(pa, pb);
            if (r < 1e-6f || r > h)
                continue;

            float3 dir = (pb - pa) / r;
            float Δ = s.restLength - r;
            float mag = k_spring * Δ * (dt * dt);
            float3 D = dir * (0.5f * mag);

            // apply equal & opposite
            newPositions[s.particleA] -= D;
            newPositions[s.particleB] += D;

            // plasticity update
            float tol = yieldRatio * s.restLength;
            if (r > s.restLength + tol)
                s.restLength += plasticity * (r - s.restLength - tol) * dt;
            else if (r < s.restLength - tol)
                s.restLength -= plasticity * (s.restLength - tol - r) * dt;

            springs[si] = s;
        }
    }
}


#endregion
