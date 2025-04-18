﻿using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using Unity.VisualScripting;
using UnityEngine;
using Unity.Burst;
using Unity.Jobs;
using Unity.Mathematics;

public class ParticleSimulator : MonoBehaviour
{
    private List<InternalParticle> internalParticles;
    private List<InternalSpring> internalSprings = new List<InternalSpring>();
    private SpatialHashGrid neighborGrid;
    public float gridCellSize = 2.0f;
    public BubbleParticleSpawner spawner;

    public Transform cameraTransform;
    public Transform mouseScreenTransform;

    private bool spawned = false;

    public Mesh particleMesh;
    public Material particleMaterial;

    public float particleRenderRadius = 0.1f;

    public float k_p = 0.004f;   // Pressure stiffness coefficient
    public float n0 = 10f;       // Rest (target) density
    public float k_p2 = 0.01f;   // Near-pressure coefficient

    public float springStiffness = 0.3f;  // k_spring, controls how strong the spring force is.
    public float plasticity = 0.3f;       // α, the plasticity constant for rest length adjustment.
    public float yieldRatio = 0.1f;       // γ, the yield threshold ratio.
    private int lastSpringUpdateFrame = -1; // to limit updates to once per frame
    public float springUpdateThreshold = 0.1f; // minimum movement needed to recheck springs

    private List<Matrix4x4> transformMatrices = new List<Matrix4x4>();

    public float gravityForce = 9.8f;
    public float particleMass;

    private float distanceToCamera;
    // Define a radius for influence.
    public float mouseDragRadius = 1.0f;
    public float mouseDragStrength = 1.0f;

    public float dragCoefficient = 0.1f; // Tweak this parameter

    private void Start()
    {
        neighborGrid = new SpatialHashGrid(gridCellSize);
        distanceToCamera = Mathf.Abs(cameraTransform.position.z - mouseScreenTransform.position.z);
    }

    private void UpdateGrid()
    {
        if(internalParticles != null)
        {
            if (neighborGrid != null)
            {
                neighborGrid.ClearGrid();
            }

            for (int i = 0; i < internalParticles.Count; i++)
            {
                neighborGrid.InsertParticle(internalParticles[i].position, i);
            }
        }
    }

    public void SpawnInternalParticles()
    {
        if (internalParticles != null) internalParticles.Clear();
        internalParticles = spawner.SpawnInternalParticles();
        UpdateGrid();

    }


    private void Update()
    {

        if(internalParticles != null)
        {
            { // rendering
                // Update your particle simulation here (or elsewhere) and update the particles list

                // Prepare instance matrices for each particle.
                int count = internalParticles.Count;
                Matrix4x4[] matrices = new Matrix4x4[count];
                for (int i = 0; i < count; i++)
                {
                    // Here we create a transformation matrix for each particle.
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
                    // 1.Prediction Phase: Apply External Forces and Predict New Positions.
                    // ---------------------------
                    // Save the current position as the old position for velocity update later.
                    for (int i = 0; i < internalParticles.Count; i++)
                    {
                        internalParticles[i].previousPosition = internalParticles[i].position;
                    }

                    // Apply external forces such as gravity to update velocity.
                    // (Here, gravityForce is assumed to be a vector; adjust if it's a scalar.)
                    for (int i = 0; i < internalParticles.Count; i++)
                    {
                        InternalParticle p = internalParticles[i];
                        // For example, adding gravity: p.velocity += gravityForce * subDeltaTime;
                        // You can also include any other external force here.
                        p.velocity += Vector3.down*gravityForce * subDeltaTime;
                        // Predict new positions by simple explicit integration:
                        p.position += p.velocity * subDeltaTime;
                    }

                    // Cache neighbor lists for all particles in this substep.
                    List<int>[] cachedNeighbors = new List<int>[internalParticles.Count];
                    for (int i = 0; i < internalParticles.Count; i++)
                    {
                        cachedNeighbors[i] = neighborGrid.GetNeighbors(internalParticles[i].position);
                    }


                    // Compute density, near-density, and pressure for each particle.
                    for (int i = 0; i < internalParticles.Count; i++)
                    {
                        ComputeDensityAndPressure(internalParticles[i], cachedNeighbors[i]);
                    }


                    // Apply displacement (position corrections) based on neighbor distances, density, etc.
                    for (int i = 0; i < internalParticles.Count; i++)
                    {
                        foreach (int neighborIndex in cachedNeighbors[i])
                        {
                            ApplyDisplacement(internalParticles[i], internalParticles[neighborIndex], subDeltaTime);
                        }
                    }

                    // Apply spring displacements:
                    for (int i = 0; i < internalSprings.Count; i++)
                    {
                        ApplySpringDisplacement(internalSprings[i], subDeltaTime);
                    }

                    // Adjust spring rest lengths:
                    for (int i = 0; i < internalSprings.Count; i++)
                    {
                        AdjustSpring(internalSprings[i], subDeltaTime);
                    }

                    // Now update each particle's velocity based on its net displacement in this substep.
                    for (int i = 0; i < internalParticles.Count; i++)
                    {
                        InternalParticle p = internalParticles[i];
                        p.velocity = (p.position - p.previousPosition) / subDeltaTime;
                    }
                }

                for (int i = 0; i < internalParticles.Count; i++)
                {
                    List<int> neighbors = neighborGrid.GetNeighbors(internalParticles[i].position);
                    foreach (int neighborIndex in neighbors)
                    {
                        TryAddSpring(i, neighborIndex);
                    }
                }


                // update the particle grid with the new positions
                UpdateGrid();

                UpdateSpringNetwork();

                ApplyMouseDragForce(Time.deltaTime);

                //ApplyAirResistance(Time.deltaTime);
            }
            
        }

        
    }

    #region Density & Pressure Computation
    private void ComputeDensityAndPressure(InternalParticle particle, List<int> neighborIndices)
    {
        // For each particle i:
        float density = 0;
        float nearDensity = 0;
        for(int i = 0; i < neighborIndices.Count;i++)
        {
            float r = Vector3.Distance(particle.position, internalParticles[neighborIndices[i]].position);
            float influenceRadius = particle.influenceRadius;
            if (r <= influenceRadius)
            {
                float q = r / influenceRadius;
                float weight = 1 - q;
                density += weight * weight;
                nearDensity += weight * weight * weight;
            }
        }
        
        particle.density = density;
        // Compute pressure:
        particle.pressure = k_p * (density - n0);
        particle.nearPressure = k_p2 * nearDensity;
    }
    #endregion

    #region Particle Displacement/ Position Update
    private void ApplyDisplacement(InternalParticle particle, InternalParticle neighbor, float subDeltaTime)
    {
        float r = Vector3.Distance(particle.position, neighbor.position);
        if (r >= particle.influenceRadius)
            return;

        float q = r / particle.influenceRadius;
        Vector3 direction = (neighbor.position - particle.position).normalized;

        // Use both particles' pressures for symmetry.
        float totalPressure = particle.pressure + neighbor.pressure;
        float totalNearPressure = particle.nearPressure + neighbor.nearPressure;
        Vector3 D = (subDeltaTime * subDeltaTime) *
                    (totalPressure * (1 - q) + totalNearPressure * Mathf.Pow(1 - q, 2)) * direction;

        // Apply equal and opposite displacements.
        particle.position -= D * 0.5f;
        neighbor.position += D * 0.5f;
    }
    #endregion

    #region Spring Forces
    private bool SpringExists(int indexA, int indexB)
    {
        return false;
        //return internalSprings.Find(s => (s.particleA == indexA && s.particleB == indexB) || (s.particleA == indexB && s.particleB == indexA)) != null;
    }

    private void TryAddSpring(int indexA, int indexB)
    {
        if (!SpringExists(indexA, indexB))
        {
            float distance = Vector3.Distance(internalParticles[indexA].position, internalParticles[indexB].position);
            // Only add the spring if the particles are sufficiently close.
            if (distance < internalParticles[indexA].influenceRadius)
            {
                internalSprings.Add(new InternalSpring(indexA, indexB, distance));
            }
        }
    }

    private void ApplySpringDisplacement(InternalSpring spring, float subDeltaTime)
    {
        InternalParticle pa = internalParticles[spring.particleA];
        InternalParticle pb = internalParticles[spring.particleB];
        Vector3 diff = pb.position - pa.position;
        float r = diff.magnitude;
        if (r >= pa.influenceRadius)
            return; // Ignore if particles are too far apart.

        Vector3 direction = diff.normalized;
        // Calculate the displacement impulse:
        Vector3 D = (subDeltaTime * subDeltaTime) * springStiffness * (1 - spring.restLength / pa.influenceRadius)
                    * (spring.restLength - r) * direction;

        // Apply equal and opposite displacements:
        pa.position -= D * 0.5f;
        pb.position += D * 0.5f;
    }

    private void AdjustSpring(InternalSpring spring, float subDeltaTime)
    {
        InternalParticle pa = internalParticles[spring.particleA];
        InternalParticle pb = internalParticles[spring.particleB];
        float r = Vector3.Distance(pa.position, pb.position);
        float d = yieldRatio * spring.restLength;  // Tolerable deformation

        if (r > spring.restLength + d)
        {
            // Spring is stretched beyond the yield; update rest length (increase it)
            spring.restLength += subDeltaTime * plasticity * (r - spring.restLength - d);
        }
        else if (r < spring.restLength - d)
        {
            // Spring is compressed beyond the yield; update rest length (decrease it)
            spring.restLength -= subDeltaTime * plasticity * (spring.restLength - d - r);
        }
    }

    private void UpdateSpringNetwork()
    {
        // Only update if we haven't already updated this frame
        if (Time.frameCount == lastSpringUpdateFrame)
            return;

        // Optionally, you might choose not to clear all springs if you want to keep persistent ones.
        // internalSprings.Clear();

        for (int i = 0; i < internalParticles.Count; i++)
        {
            // Only update springs for particle i if it moved more than the threshold.
            if (Vector3.Distance(internalParticles[i].position, internalParticles[i].lastSpringUpdatePosition) > springUpdateThreshold)
            {
                // Update its last checked position.
                internalParticles[i].lastSpringUpdatePosition = internalParticles[i].position;

                // Get neighbors for particle i.
                List<int> neighbors = neighborGrid.GetNeighbors(internalParticles[i].position);
                foreach (int neighborIndex in neighbors)
                {
                    // Only add each pair once.
                    if (neighborIndex > i)
                    {
                        TryAddSpring(i, neighborIndex);
                    }
                }
            }
        }

        lastSpringUpdateFrame = Time.frameCount;
    }

    #endregion

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
}


