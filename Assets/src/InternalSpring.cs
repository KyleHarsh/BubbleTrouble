using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class InternalSpring 
{
    public int particleA;
    public int particleB;
    public float restLength;

    public InternalSpring(int a, int b, float restLength)
    {
        this.particleA = a;
        this.particleB = b;
        this.restLength = restLength;
    }
}
