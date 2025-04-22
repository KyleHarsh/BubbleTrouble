using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

public class InternalParticle
{
    public float influenceRadius = 1.0f;
    public Vector3 normal;
    public Vector3 position;
    public Vector3 previousPosition;
    public Vector3 velocity;
    public float density;
    public float nearDensity;
    public float pressure;
    public float nearPressure;
    public Vector3 lastSpringUpdatePosition = Vector3.zero;


    public InternalParticle(float density, Vector3 position, Vector3 velocity, float influence)
    {
        this.position = position;
        this.velocity = velocity;
        this.density = density;
        influenceRadius = influence;
    }
}

public struct NeighborPair
{
    public int i, j;
    public NeighborPair(int i, int j)
    {
        this.i = i;
        this.j = j;
    }
}

