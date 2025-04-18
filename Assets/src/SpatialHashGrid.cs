using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Collections;
using Unity.Mathematics;
using Unity.Jobs;
using Unity.Burst;

public class SpatialHashGrid 
{
    private Dictionary<Vector3Int, List<int>> grid = new Dictionary<Vector3Int, List<int>>();
    private float cellSize;

    public SpatialHashGrid(float cellSize)
    {
        this.cellSize = cellSize;
    }

    public Vector3Int GetCell(Vector3 position)
    {
        return new Vector3Int(
            Mathf.FloorToInt(position.x / cellSize),
            Mathf.FloorToInt(position.y / cellSize),
            Mathf.FloorToInt(position.z / cellSize)
        );
    }

    public void ClearGrid()
    {
        grid.Clear();
    }

    public void InsertParticle(Vector3 position, int index)
    {
        Vector3Int cell = GetCell(position);
        if (!grid.ContainsKey(cell))
        {
            grid[cell] = new List<int>();
        }
        grid[cell].Add(index);
    }

    // This method queries all particles in neighboring cells
    public List<int> GetNeighbors(Vector3 position)
    {
        List<int> neighbors = new List<int>();
        Vector3Int baseCell = GetCell(position);

        for (int x = -1; x <= 1; x++)
        {
            for (int y = -1; y <= 1; y++)
            {
                for (int z = -1; z <= 1; z++)
                {
                    Vector3Int neighborCell = baseCell + new Vector3Int(x, y, z);
                    if (grid.ContainsKey(neighborCell))
                    {
                        neighbors.AddRange(grid[neighborCell]);
                    }
                }
            }
        }
        return neighbors;
    }
}

public struct NativeHashGrid : System.IDisposable
{
    public NativeParallelMultiHashMap<int, int> map;
    private float cellSize;

    public NativeHashGrid(int capacity, float cellSize, Allocator alloc)
    {
        map = new NativeParallelMultiHashMap<int, int>(capacity, alloc);
        this.cellSize = cellSize;
    }

    public void Clear() => map.Clear();

    private static int HashCell(int x, int y, int z)
    {
        unchecked
        {
            int h = x * 73856093;
            h ^= y * 19349663;
            h ^= z * 83492791;
            return h;
        }
    }

    public void Insert(float3 pos, int index)
    {
        int x = (int)math.floor(pos.x / cellSize);
        int y = (int)math.floor(pos.y / cellSize);
        int z = (int)math.floor(pos.z / cellSize);
        int key = HashCell(x, y, z);
        map.Add(key, index);
    }

    /// <summary>
    /// Job‐safe neighbors query. 
    /// Callers must supply a NativeList<int> with a Temp or TempJob allocator.
    /// </summary>
    public void GetNeighbors(float3 pos, ref NativeList<int> results)
    {
        results.Clear();
        int cx = (int)math.floor(pos.x / cellSize);
        int cy = (int)math.floor(pos.y / cellSize);
        int cz = (int)math.floor(pos.z / cellSize);

        // loop 3×3×3 cells
        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                {
                    int key = HashCell(cx + dx, cy + dy, cz + dz);
                    if (map.TryGetFirstValue(key, out int neighbor, out var iter))
                    {
                        do
                        {
                            results.Add(neighbor);
                        } while (map.TryGetNextValue(out neighbor, ref iter));
                    }
                }
    }

    public void Dispose() => map.Dispose();
}
