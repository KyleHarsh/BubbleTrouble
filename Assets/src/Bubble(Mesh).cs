using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Unity.VisualScripting;
using UnityEditor.PackageManager;
using UnityEngine;

public class Bubble_m : MonoBehaviour
{
    // Start is called before the first frame update

    private MeshFilter meshF;
    private List<SurfaceParticle> surfaceParticles;
    private int numTries;

    public int numSimsPerFrame = 3;
    public float springConst = 0.1f;
    public float damping = 0.01f;

    public Camera cam;
    public float clickRadius = 0.7f;
    public float clickStrength = 10.5f;

    private bool flag = false;

    private Vector3 initialVelocity = Vector3.zero;
    private Vector3 initialAcceleration = Vector3.zero;
    void Start()
    {
        meshF = GetComponent<MeshFilter>();
        //LazySquirrelLabs.SphereGenerator.Generators.IcosphereGenerator gen = new LazySquirrelLabs.SphereGenerator.Generators.IcosphereGenerator(1.0f, 1);
        //this.meshF.mesh = gen.Generate();
        SphereGenerator sphereGenerator = new SphereGenerator();
        this.meshF.mesh = sphereGenerator.Create(3);
        this.surfaceParticles = new List<SurfaceParticle>();
        removeDuplicateVerticies();
        meshToBubble();
        numTries = this.meshF.mesh.triangles.Length;
    }
    void Update()
    {
        
        for (int i = 0; i < numSimsPerFrame; ++i)
        {
            runSim();
        }


        bubbleToMesh();
        CheckMouseClick();

    }

    #region mouse
    void CheckMouseClick()
    {
        if (Input.GetMouseButton(0))
        {  // If the mouse button is pressed
            Debug.Log("Mouse Clicked");
            Ray ray = cam.ScreenPointToRay(Input.mousePosition);
            //for each particle
            for (int i = 0; i < surfaceParticles.Count; ++i)
            {
                //if the particle is close to the ray cast by the mouse, give it a push in that direction
                float dist = Vector3.Cross(ray.direction, surfaceParticles[i].getPosition() - ray.origin).magnitude;
                if (dist <= clickRadius)
                {
                    surfaceParticles[i].applyForce(ray.direction*clickStrength);
                }
            }
        }
    }

    #endregion
    #region sim
    private Vector3 forceOnParticle(int i)
    {
        List<int> neighbors = this.surfaceParticles[i].getNeigbors();
        List<float> springDists = this.surfaceParticles[i].getDists();
        Vector3 force = Vector3.zero;
        for (int j = 0; j < neighbors.Count; ++j)
        {
            Vector3 diff = surfaceParticles[neighbors[j]].getPosition() - surfaceParticles[i].getPosition();
            Vector3 dir = Vector3.Normalize(diff);
            float dist = diff.magnitude;
            force += -1 * springConst * (springDists[j] -dist) * 0.5f * dir;
        }
        return force;
    }
    private void runSim()
    {
        for (int i = 0; i < surfaceParticles.Count; i++)
        {
            surfaceParticles[i].applyForce(forceOnParticle(i));
            surfaceParticles[i].update(Time.deltaTime/numSimsPerFrame, damping);
        }
    }
    #endregion
    #region dupicates
    
    private void removeDuplicateVerticies()
    {
        Dictionary<Vector3, int> duplicates = new Dictionary<Vector3, int>();
        Vector3[] verts = meshF.mesh.vertices;
        for (int i = 0; i < verts.Length; i++)
        {
            if (!duplicates.ContainsKey(verts[i]))
            {
                duplicates[verts[i]] = -1;

            } else
            {
                //Debug.Log(string.Format("DUPLICATE FOUND: {0}", verts[i].ToString()));
            }
        }

        //makes the duplicates int it's new place
        Vector3[] unDupedVerts = duplicates.Keys.ToArray();
        for (int i = 0;i < unDupedVerts.Length; i++)
        {
            duplicates[unDupedVerts[i]] = i;
        }

        int[] tries = this.meshF.mesh.triangles;
        int[] unDupedTries = new int[tries.Length];
        for (int i = 0; i < tries.Length; i++)
        {
            unDupedTries[i] = duplicates[verts[tries[i]]];
        }
        this.meshF.mesh.triangles = unDupedTries;
        this.meshF.mesh.vertices = unDupedVerts;
        Debug.Log(string.Join(", ", unDupedVerts));
        Debug.Log(string.Format("Verticies removed: {0}", verts.Length - unDupedVerts.Length));
        //Debug.Log(string.Format("Verticies length: {0}", unDupedVerts.Length));
        //Debug.Log(string.Format("Max index in tries: {0}", Enumerable.Max(unDupedTries)));
        //Debug.Log(string.Join(", ", unDupedTries));
    }
    #endregion
    #region construction

    void meshToBubble()
    {
        if (surfaceParticles.Count != 0)
        {
            surfaceParticles.Clear();
        }
        int[] tries = meshF.mesh.triangles;
        Vector3[] verts = meshF.mesh.vertices;

        //I am assuming we are using object coordinates
        for (int i  = 0; i < verts.Length; ++i)
        {
            surfaceParticles.Add(new SurfaceParticle(verts[i], i));
        }
        //do neigbors
        for (int i = 0; i < tries.Length; i += 3)
        {
            float dist01 = Math.Abs((surfaceParticles[tries[i]].getPosition() - surfaceParticles[tries[i+1]].getPosition()).magnitude);
            float dist12 = Math.Abs((surfaceParticles[tries[i+1]].getPosition() - surfaceParticles[tries[i + 2]].getPosition()).magnitude);
            float dist02 = Math.Abs((surfaceParticles[tries[i]].getPosition() - surfaceParticles[tries[i + 2]].getPosition()).magnitude);

            //only one of the verticies knows about the triangle so that we don't have duplicate triangles on mesh construction
            surfaceParticles[tries[i]].addNeighbor(tries[i + 1], dist01, tries[i + 2]);
            surfaceParticles[tries[i]].addNeighbor(tries[i + 2], dist02);

            surfaceParticles[tries[i + 1]].addNeighbor(tries[i], dist01);
            surfaceParticles[tries[i + 1]].addNeighbor(tries[i + 2], dist12);

            surfaceParticles[tries[i + 2]].addNeighbor(tries[i], dist02);
            surfaceParticles[tries[i + 2]].addNeighbor(tries[i + 1], dist12);
        }
    }
    void bubbleToMesh()
    {
        Vector3[] verts = new Vector3[this.surfaceParticles.Count];
        for (int i = 0; i < this.surfaceParticles.Count; ++i)
        {
            verts[i] = (this.surfaceParticles[i].getPosition());
        }

        int[] tries = new int[this.numTries];

        int t = 0;

        for (int i = 0; i < this.surfaceParticles.Count; ++i)
        {
            List<(int, int)> iTr = this.surfaceParticles[i].getTries();
            for (int j = 0; j < iTr.Count; ++j)
            {
                if (t > this.numTries - 3)
                {
                    print("MESH CONSTRUCTION ERROR: too many trinagles");
                }
                tries[t] = i;
                tries[t + 1] = iTr[j].Item1;
                tries[t + 2] = iTr[j].Item2;
                t+=3;
            }
        }
        this.meshF.mesh.vertices = verts;
        this.meshF.mesh.triangles = tries;
    }
    #endregion
}
