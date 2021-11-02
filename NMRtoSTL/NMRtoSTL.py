from io import FileIO
import os
from typing import Tuple
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as triang
import numpy as np
import re
from stl import mesh
import nmrglue as ng
from stl.stl import ASCII
import importNMR


def create_mesh(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, verbose=False
) -> mesh.Mesh:
    if verbose:
        print("Triangulating")
    tri = triang.Triangulation(
        y,
        x,
    )
    data = np.zeros(len(tri.triangles), dtype=mesh.Mesh.dtype)
    print("Constructing Mesh")
    NMR_mesh = mesh.Mesh(data, remove_empty_areas=False)
    NMR_mesh.x = x[tri.triangles]
    NMR_mesh.y = y[tri.triangles]
    NMR_mesh.z = z[tri.triangles]
    return NMR_mesh


def create_base(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    thickness=1,
):
    """
    Creates a base for the 3d model, uses code from the numpy-STL documentation
    to create the base:
    https://numpy-stl.readthedocs.io/en/latest/usage.html#creating-mesh-objects-from-a-list-of-vertices-and-faces
    """
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    zmin = z.min()

    # Create the vertices matrix for the base
    # the min and max of the x and y make the base as long and wide as the data
    # the zmin starts the base at the lowest point in the data and makes it
    # as thick as the thickness.
    vertices = np.array(
        [
            [xmin, ymin, zmin - thickness],
            [xmax, ymin, zmin - thickness],
            [xmax, ymax, zmin - thickness],
            [xmin, ymax, zmin - thickness],
            [xmin, ymin, zmin],
            [xmax, ymin, zmin],
            [xmax, ymax, zmin],
            [xmin, ymax, zmin],
        ]
    )
    # Define the 12 triangles composing the cube
    faces = np.array(
        [
            [0, 3, 1],
            [1, 3, 2],
            [0, 4, 7],
            [0, 7, 3],
            [4, 5, 6],
            [4, 6, 7],
            [5, 1, 2],
            [5, 2, 6],
            [2, 3, 6],
            [3, 7, 6],
            [0, 1, 5],
            [0, 5, 4],
        ]
    )

    # Create the mesh
    base = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            base.vectors[i][j] = vertices[f[j], :]

    return base


def main(filename):
    spectrum = importNMR.importNMR(filename)
    spectrum.read_file(verbose=True)
    spectrum.make_scales(verbose=True, 
                         f1_min=1.0, f1_max=4.5, f2_min=1.0, f2_max=4.5)
    x, y, z = spectrum.process()
    NMR_mesh = create_mesh(x, y, z, verbose=True)
    _, file = os.path.split(filename)
    file, _ = os.path.splitext(file)
    filename = f"./{file}.stl"
    print(f"Saving file as: '{filename[2:]}'")
    NMR_mesh.save(filename)


if __name__ == "__main__":
    import sys

    try:
        filename = sys.argv[1]
    except:
        filename = r"..\Example_data\Bruker_COSY\pdata\1"
    main(filename)
