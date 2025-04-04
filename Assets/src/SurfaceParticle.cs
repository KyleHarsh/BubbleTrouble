using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;

public class SurfaceParticle
{

    public SurfaceParticle(UnityEngine.Vector3 p, UnityEngine.Vector3 v, UnityEngine.Vector3 a, float m = 0.001f)
    {
        this.position = p;
        this.velocity = v;
        this.acceleration = a;
        this.mass = m;
        this.neigbors = new List<SurfaceParticle>();
    }

    public void addNeighbor(SurfaceParticle n)
    {
        this.neigbors.Add(n);
    }

    public UnityEngine.Vector3 getPosition() { return this.position; }
    public void externalForce(UnityEngine.Vector3 f)
    {
        acceleration += f * this.mass;
    }

    // Update is called once per frame
    void Update()
    {
        position += velocity * Time.deltaTime;
        velocity += acceleration * Time.deltaTime;

        // surface tension kinda
        // particles are assumed to be 0.1 unit away from each other by default
        for (int i = 0; i < neigbors.Count; i++)
        {
            UnityEngine.Vector3 direction = (neigbors[i].position - this.position).normalized;
            acceleration -= k * ((neigbors[i].position - this.position).magnitude - partDist) * direction;
        }
        acceleration -= damping * Time.deltaTime * velocity;
    }

    private static readonly float k = 0.1f;
    private static readonly float partDist = 0.1f;
    private static readonly float damping = 0.1f; //this is the ammount of damping/sec

    //representation:
    float mass;
    UnityEngine.Vector3 position;
    UnityEngine.Vector3 velocity;
    UnityEngine.Vector3 acceleration;

    List<SurfaceParticle> neigbors;
}
