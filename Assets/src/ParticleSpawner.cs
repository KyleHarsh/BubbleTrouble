using System;
using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.Scripting;

public class BubbleParticleSpawner : MonoBehaviour
{
    [SerializeField] private Transform spawnCenter;
    public BubbleParticleType typeToSpawn;
    public int numParticlesToSpawn = 10;
    public float minParticleDensity = 4;
    public float bubbleRadius = 3.0f;

    List<InternalParticle> particlesSpawned = new List<InternalParticle>();

    public static event Action OnParticlesSpawned;

    private void Start()
    {
        if (spawnCenter == null) spawnCenter = transform;
    }
    public List<InternalParticle> SpawnInternalParticles()
    {
        if (particlesSpawned != null) particlesSpawned.Clear();

        particlesSpawned = new List<InternalParticle> ();

        // Calculate the volume of the sphere
        float sphereVolume = (4f / 3f) * Mathf.PI * Mathf.Pow(bubbleRadius, 3);

        // Calculate the average volume per point and then the spacing
        float avgVolumePerPoint = sphereVolume / numParticlesToSpawn;
        float spacing = Mathf.Pow(avgVolumePerPoint, 1f / 3f);

        List<Vector3> particlePositions = PoissonDiskSampler3D.GeneratePoints(spacing, bubbleRadius, numParticlesToSpawn);

        for(int i = 0; i < particlePositions.Count; i++)
        {
            InternalParticle particle = new InternalParticle(minParticleDensity, particlePositions[i]+transform.position, Vector3.zero, minParticleDensity);
            particlesSpawned.Add(particle);

        }

        return particlesSpawned;

    }

    private void OnDrawGizmos()
    {
        Gizmos.color = Color.blue.WithAlpha(0.5f);
        Gizmos.DrawWireSphere(spawnCenter.position, bubbleRadius);
    }
}

public enum BubbleParticleType 
{ 
    Surface,
    Internal
}
