#%%
from io import FileIO
from typing import Union
import matplotlib.pyplot as plt
import matplotlib.tri as triang
import numpy as np
import re
from stl import mesh
ROWS_RE = re.compile(r"NROWS = (\d+)")
COLS_RE = re.compile(r"NCOLS = (\d+)")
XY_RE = re.compile(r"F[12][\w]{4,5} = ([-\d.]+)")

def get_consts(file: FileIO):
    # read the head of the file to find the number of rows, columns and 
    # the limits of the X and Y axes
    header = file.read(600)
    file.seek(0)
    #find rows and columns
    rows = int(re.findall(ROWS_RE, header)[0])
    cols = int(re.findall(COLS_RE, header)[0])
    #find axis limits
    xmax, xmin, ymax, ymin = [float(value) for value in re.findall(XY_RE, header)]
    return rows, cols, (xmax, xmin, ymax, ymin)

def read_file(filename, verbose = False):
    if verbose:
        print("Reading Data")
    # reads data and constants from the text file 
    data = []
    with open(filename, "r") as f:
        NROWS, NCOLS, axes = get_consts(f)
        for line in f:
            if line.startswith("#"):
                continue
            else:
                row = [float(line)] + [float(f.readline().strip()) for _ in range(NCOLS-1)]
                data.append(row)
    return data, axes, (NROWS, NCOLS)

def main(filename):
    data, axes, shape = read_file(filename, verbose=True)
    xmax, xmin, ymax, ymin = axes
    NROWS, NCOLS = shape
    x = np.linspace(xmax, xmin, num=NROWS)
    y = np.linspace(ymax, ymin, num=NCOLS)

    arr =np.array(data)
    X, Y = np.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
    z = arr.flatten()
    # normalise z data
    z = (z/z.max()) * 1
    x = X.flatten()
    y = Y.flatten()
    print("performing Triangulation")
    tri=triang.Triangulation(y,x)
    data = np.zeros(len(tri.triangles), dtype=mesh.Mesh.dtype)
    print("Constructing Mesh")
    mobius_mesh = mesh.Mesh(data, remove_empty_areas=False)
    mobius_mesh.x = x[tri.triangles]
    mobius_mesh.y[:] = y[tri.triangles]
    mobius_mesh.z[:] = z[tri.triangles]
    mobius_mesh.save(f'{filename}.stl')

if __name__ == "__main__":
    import sys
    try:
        filename = sys.argv[1]
    except:
        filename = "NMR2.txt"
    main(filename)