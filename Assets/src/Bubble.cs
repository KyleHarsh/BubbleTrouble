using System.Collections;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using UnityEngine;


public class Bubble : MonoBehaviour
{
    void spawnParticles(int num, UnityEngine.Vector3 center)
    {
        for (int i = 0; i < num; i++)
        {
            UnityEngine.Vector3 pos = Random.onUnitSphere*this.radius;
            SurfaceParticle particle = new SurfaceParticle(center + pos, UnityEngine.Vector3.zero, UnityEngine.Vector3.zero);
            this.surfaceParticles.Add(particle);
        }
        for (int i = 0; i < surfaceParticles.Count; ++i)
        {
            for(int j = 0; j < surfaceParticles.Count; ++j)
            {
                if (i == j) break;
                if ((surfaceParticles[i].getPosition() - surfaceParticles[j].getPosition()).magnitude <= neighborDist)
                {
                    surfaceParticles[i].addNeighbor(surfaceParticles[j]);
                    surfaceParticles[j].addNeighbor(surfaceParticles[i]);
                }
            }
        }
    }

    // Start is called before the first frame update
    void Start()
    {
        this.spawnParticles(100, UnityEngine.Vector3.zero);
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    UnityEngine.Vector3 calcCenterOfMass()
    {
        //not yet implemented
        return UnityEngine.Vector3.zero;
    }

    private static readonly float neighborDist = 0.1f;

    List<SurfaceParticle> surfaceParticles;
    float radius;
}
