using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;

public class Bubble_m : MonoBehaviour
{
    // Start is called before the first frame update

    private MeshFilter meshF;
    private List<SurfaceParticle> surfaceParticles;

    private UnityEngine.Vector3 initialVelocity = UnityEngine.Vector3.zero;
    private UnityEngine.Vector3 initialAcceleration = UnityEngine.Vector3.zero;
    void Start()
    {
        meshF = GetComponent<MeshFilter>();
        LazySquirrelLabs.SphereGenerator.Generators.IcosphereGenerator gen = new LazySquirrelLabs.SphereGenerator.Generators.IcosphereGenerator(1.0f, 3);
        this.meshF.mesh = gen.Generate();
        meshToBubble();
    }

    void meshToBubble()
    {
        if (surfaceParticles.Count != 0)
        {
            surfaceParticles.Clear();
        }
        int[] tries = meshF.mesh.triangles;
        UnityEngine.Vector3[] verts = meshF.mesh.vertices;

        Dictionary<int, int> partIndicies = new Dictionary<int, int>();
        //I am assuming we are using world coordinates
        for (int i = 0; i < tries.Length; i += 3)
        {
            UnityEngine.Vector3 pos1 = transform.TransformPoint(verts[tries[i]]);
            UnityEngine.Vector3 pos2 = transform.TransformPoint(verts[tries[i+1]]);
            UnityEngine.Vector3 pos3 = transform.TransformPoint(verts[tries[i+2]]);

            int part1;
            int part2;
            int part3;

            if (!partIndicies.ContainsKey(tries[i]))
            {
                surfaceParticles.Add(new SurfaceParticle(pos1, initialVelocity, initialAcceleration));
                partIndicies.Add(tries[i], surfaceParticles.Count - 1);
            }
            part1 = partIndicies[tries[i]];
            if (!partIndicies.ContainsKey(tries[i+1]))
            {
                surfaceParticles.Add(new SurfaceParticle(pos2, initialVelocity, initialAcceleration));
                partIndicies.Add(tries[i+1], surfaceParticles.Count - 1);
            }
            part2 = partIndicies[tries[i+1]];
            if (!partIndicies.ContainsKey(tries[i+2]))
            {
                surfaceParticles.Add(new SurfaceParticle(pos3, initialVelocity, initialAcceleration));
                partIndicies.Add(tries[i+2], surfaceParticles.Count - 1);
            }
            part3 = partIndicies[tries[i+2]];

            surfaceParticles[part1].addNeighbor(surfaceParticles[part2]);
            surfaceParticles[part1].addNeighbor(surfaceParticles[part3]);

            surfaceParticles[part2].addNeighbor(surfaceParticles[part1]);
            surfaceParticles[part2].addNeighbor(surfaceParticles[part3]);

            surfaceParticles[part3].addNeighbor(surfaceParticles[part1]);
            surfaceParticles[part3].addNeighbor(surfaceParticles[part2]);
        }
    }
    void bubbleToMesh()
    {

    }
}
