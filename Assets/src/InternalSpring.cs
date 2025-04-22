using System.Collections;
using System.Collections.Generic;
using UnityEngine;
public struct InternalSpring : IEnumerable<InternalSpring>
{
    public int particleA;
    public int particleB;
    public float restLength;

    public InternalSpring(int a, int b, float restLength)
    {
        particleA = a;
        particleB = b;
        this.restLength = restLength;
    }

    public IEnumerator<InternalSpring> GetEnumerator()
    {
        throw new System.NotImplementedException();
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        throw new System.NotImplementedException();
    }
}

public class InternalSpringClass
{
    public int particleA;
    public int particleB;
    public float restLength;

    public InternalSpringClass(int a, int b, float restLength)
    {
        particleA = a;
        particleB = b;
        this.restLength = restLength;
    }
}
