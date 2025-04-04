using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PoissonDiskSampler3D
{
    // Returns a list of points uniformly distributed in 3D, 
    // with a minimum distance 'minDist' between them, within a sphere of radius 'domainRadius'.
    public static List<Vector3> GeneratePoints(float minDist, float domainRadius, int k = 30)
    {
        float cellSize = minDist / Mathf.Sqrt(3);
        // Compute grid dimensions for a cube that encloses the sphere.
        int gridSize = Mathf.CeilToInt((2 * domainRadius) / cellSize);
        // 3D grid represented as a dictionary mapping grid indices to sample points.
        Dictionary<Vector3Int, Vector3> grid = new Dictionary<Vector3Int, Vector3>();

        List<Vector3> samples = new List<Vector3>();
        List<Vector3> active = new List<Vector3>();

        // Helper: Convert a position to grid coordinate.
        Vector3Int GetGridPos(Vector3 pos)
        {
            // Shift the coordinate system so that the domain spans from (0,0,0) to (2*domainRadius)
            int x = Mathf.FloorToInt((pos.x + domainRadius) / cellSize);
            int y = Mathf.FloorToInt((pos.y + domainRadius) / cellSize);
            int z = Mathf.FloorToInt((pos.z + domainRadius) / cellSize);
            return new Vector3Int(x, y, z);
        }

        // Check if a candidate is far enough from existing samples.
        bool IsFarEnough(Vector3 candidate)
        {
            Vector3Int gridPos = GetGridPos(candidate);
            // Check neighboring cells in a 3x3x3 block.
            for (int x = -1; x <= 1; x++)
            {
                for (int y = -1; y <= 1; y++)
                {
                    for (int z = -1; z <= 1; z++)
                    {
                        Vector3Int neighborPos = new Vector3Int(gridPos.x + x, gridPos.y + y, gridPos.z + z);
                        if (grid.TryGetValue(neighborPos, out Vector3 neighborSample))
                        {
                            if (Vector3.Distance(candidate, neighborSample) < minDist)
                                return false;
                        }
                    }
                }
            }
            return true;
        }

        // Check if candidate is inside the spherical domain.
        bool IsInDomain(Vector3 candidate)
        {
            return candidate.sqrMagnitude <= domainRadius * domainRadius;
        }

        // Generate a random point in the spherical shell between [minDist, 2*minDist] around a given point.
        Vector3 GenerateRandomPointAround(Vector3 point)
        {
            // Random direction using a uniform point on the unit sphere.
            Vector3 dir = Random.onUnitSphere;
            // Random distance between minDist and 2*minDist.
            float distance = Random.Range(minDist, 2 * minDist);
            return point + dir * distance;
        }

        // Start with an initial random point within the domain.
        Vector3 initialPoint = Random.insideUnitSphere * domainRadius;
        samples.Add(initialPoint);
        active.Add(initialPoint);
        grid[GetGridPos(initialPoint)] = initialPoint;

        // Main loop: Process the active list.
        while (active.Count > 0)
        {
            int index = Random.Range(0, active.Count);
            Vector3 point = active[index];
            bool foundCandidate = false;
            for (int i = 0; i < k; i++)
            {
                Vector3 candidate = GenerateRandomPointAround(point);
                if (IsInDomain(candidate) && IsFarEnough(candidate))
                {
                    samples.Add(candidate);
                    active.Add(candidate);
                    grid[GetGridPos(candidate)] = candidate;
                    foundCandidate = true;
                    break;
                }
            }
            if (!foundCandidate)
            {
                active.RemoveAt(index);
            }
        }

        return samples;
    }
}
