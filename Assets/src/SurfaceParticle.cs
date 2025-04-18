using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

public class SurfaceParticle
{

    public SurfaceParticle(Vector3 p, Vector3 v, Vector3 a, int id_, float m = 0.001f)
    {
        this.position = p;
        this.velocity = v;
        this.acceleration = a;
        this.mass = m;
        this.neigbors = new List<(int, int)>();
        this.distances = new List<float> ();
        this.id = id_;
    }
    public SurfaceParticle(Vector3 p, int id_, float m = 0.001f)
    {
        this.position = p;
        this.velocity = Vector3.zero;
        this.acceleration = Vector3.zero;
        this.mass = m;
        this.neigbors = new List<(int, int)>();
        this.distances = new List<float>();
        this.id = id_;
    }

    public void addNeighbor(int a, float d, int b = -1)
    {
        this.neigbors.Add((a, b));
        this.distances.Add(d);
    }

    public Vector3 getPosition() { return this.position; }
    public List<int> getNeigbors()
    {
        List<int> ret = new List<int>();
        for (int i = 0; i < this.neigbors.Count; ++i)
        {
            ret.Add(this.neigbors[i].Item1);
        }

        return ret;
    }

    public List<float> getDists() { return this.distances; }
    public List<(int, int)> getTries()
    {
        List<(int, int)> ret = new List<(int, int)>();
        for (int i = 0; i < this.neigbors.Count; ++i)
        {
            if (this.neigbors[i].Item2 != -1)
            {
                ret.Add(this.neigbors[i]);
            }
        }

        return ret;
    }

    public void applyForce(Vector3 f)
    {
        acceleration += f * this.mass;
    }

    public void update(float dt, float damping)
    {
        this.position += this.velocity * dt;
        this.velocity += this.acceleration * dt;
        this.acceleration -= this.acceleration * dt * damping;
    }

    //representation:
    float mass;
    Vector3 position;
    Vector3 velocity;
    Vector3 acceleration;

    //ID is the index in the array it is stored in
    int id;

    //id of the neighbor and optional id of third vertex in triangle
    List<(int, int)> neigbors;

    List<float> distances;
}
