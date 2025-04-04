using System.Collections;
using System.Collections.Generic;
using UnityEngine;

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
