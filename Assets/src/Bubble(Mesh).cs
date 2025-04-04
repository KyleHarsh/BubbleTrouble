using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Bubble_m : MonoBehaviour
{
    // Start is called before the first frame update

    private MeshFilter meshF;
    void Start()
    {
        meshF = GetComponent<MeshFilter>();
        LazySquirrelLabs.SphereGenerator.Generators.IcosphereGenerator gen = new LazySquirrelLabs.SphereGenerator.Generators.IcosphereGenerator(1.0f, 3);
        this.meshF.mesh = gen.Generate();

        // create new colors array where the colors will be created.
        Color[] colors = new Color[meshF.mesh.vertices.Length];

        for (int i = 0; i < meshF.mesh.vertices.Length; i++)
            //colors[i] = Color.Lerp(Color.red, Color.blue, meshF.mesh.vertices[i].y);
            colors[i] = Color.cyan;

        // assign the array of colors to the Mesh.
        meshF.mesh.colors = colors;
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
