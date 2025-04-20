using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;

public class SphereGenerator
{
    private struct TriangleIndices
    {
        public int v1;
        public int v2;
        public int v3;

        public TriangleIndices(int v1, int v2, int v3)
        {
            this.v1 = v1;
            this.v2 = v2;
            this.v3 = v3;
        }
    }

    private int index;
    private Dictionary<Int64, int> middlePointIndexCache;
    private List<Vector3> positions;

    // add vertex to mesh, fix position to be on unit sphere, return index
    private int addVertex(Vector3 p, float r)
    {
        float length = (float)Math.Sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        this.positions.Add(new Vector3(p.x*r / length, p.y*r / length, p.z*r / length));
        return index++;
    }

    // return index of point in the middle of p1 and p2
    private int getMiddlePoint(int p1, int p2, float r)
    {
        // first check if we have it already
        bool firstIsSmaller = p1 < p2;
        Int64 smallerIndex = firstIsSmaller ? p1 : p2;
        Int64 greaterIndex = firstIsSmaller ? p2 : p1;
        Int64 key = (smallerIndex << 32) + greaterIndex;

        int ret;
        if (this.middlePointIndexCache.TryGetValue(key, out ret))
        {
            return ret;
        }

        // not in cache, calculate it
        Vector3 point1 = this.positions[p1];
        Vector3 point2 = this.positions[p2];
        Vector3 middle = new Vector3(
            (float)((point1.x + point2.x) / 2.0),
            (float)((point1.y + point2.y) / 2.0),
            (float)((point1.z + point2.z) / 2.0));

        // add vertex makes sure point is on unit sphere
        int i = addVertex(middle, r);

        // store it, return index
        this.middlePointIndexCache.Add(key, i);
        return i;
    }

    public Mesh Create(float radius, int recursionLevel)
    {
        Mesh mesh = new Mesh();
        this.middlePointIndexCache = new Dictionary<long, int>();
        this.index = 0;

        this.positions = new List<Vector3>();

        // create 12 vertices of a icosahedron
        float t = (float)((1.0 + Math.Sqrt(5.0)) / 2.0);

        addVertex(new Vector3(-1, t, 0), radius);
        addVertex(new Vector3(1, t, 0), radius);
        addVertex(new Vector3(-1, -t, 0), radius);
        addVertex(new Vector3(1, -t, 0), radius);

        addVertex(new Vector3(0, -1, t), radius);
        addVertex(new Vector3(0, 1, t), radius);
        addVertex(new Vector3(0, -1, -t), radius);
        addVertex(new Vector3(0, 1, -t), radius);

        addVertex(new Vector3(t, 0, -1), radius);
        addVertex(new Vector3(t, 0, 1), radius);
        addVertex(new Vector3(-t, 0, -1), radius);
        addVertex(new Vector3(-t, 0, 1), radius);


        // create 20 triangles of the icosahedron
        var faces = new List<TriangleIndices>();

        // 5 faces around point 0
        faces.Add(new TriangleIndices(0, 11, 5));
        faces.Add(new TriangleIndices(0, 5, 1));
        faces.Add(new TriangleIndices(0, 1, 7));
        faces.Add(new TriangleIndices(0, 7, 10));
        faces.Add(new TriangleIndices(0, 10, 11));

        // 5 adjacent faces 
        faces.Add(new TriangleIndices(1, 5, 9));
        faces.Add(new TriangleIndices(5, 11, 4));
        faces.Add(new TriangleIndices(11, 10, 2));
        faces.Add(new TriangleIndices(10, 7, 6));
        faces.Add(new TriangleIndices(7, 1, 8));

        // 5 faces around point 3
        faces.Add(new TriangleIndices(3, 9, 4));
        faces.Add(new TriangleIndices(3, 4, 2));
        faces.Add(new TriangleIndices(3, 2, 6));
        faces.Add(new TriangleIndices(3, 6, 8));
        faces.Add(new TriangleIndices(3, 8, 9));

        // 5 adjacent faces 
        faces.Add(new TriangleIndices(4, 9, 5));
        faces.Add(new TriangleIndices(2, 4, 11));
        faces.Add(new TriangleIndices(6, 2, 10));
        faces.Add(new TriangleIndices(8, 6, 7));
        faces.Add(new TriangleIndices(9, 8, 1));


        // refine triangles
        for (int i = 0; i < recursionLevel; i++)
        {
            var faces2 = new List<TriangleIndices>();
            foreach (var tri in faces)
            {
                // replace triangle by 4 triangles
                int a = getMiddlePoint(tri.v1, tri.v2, radius);
                int b = getMiddlePoint(tri.v2, tri.v3, radius);
                int c = getMiddlePoint(tri.v3, tri.v1, radius);

                faces2.Add(new TriangleIndices(tri.v1, a, c));
                faces2.Add(new TriangleIndices(tri.v2, b, a));
                faces2.Add(new TriangleIndices(tri.v3, c, b));
                faces2.Add(new TriangleIndices(a, b, c));
            }
            faces = faces2;
        }

        List<int> triangles = new List<int>();

        // done, now add triangles to mesh
        foreach (var tri in faces)
        {
            triangles.Add(tri.v1);
            triangles.Add(tri.v2);
            triangles.Add(tri.v3);
        }

        mesh.SetVertices(this.positions);
        mesh.SetIndices(triangles, MeshTopology.Triangles, 0);
        mesh.RecalculateBounds();
        mesh.RecalculateNormals();
        mesh.RecalculateTangents();

        return mesh;
    }
}
