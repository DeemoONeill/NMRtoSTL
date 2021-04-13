#%%
from io import FileIO
import os
from typing import Union
import matplotlib.pyplot as plt
import matplotlib.tri as triang
import numpy as np
import re
from stl import mesh


class Bruker2D:

    ROWS_RE = re.compile(r"NROWS = (\d+)")
    COLS_RE = re.compile(r"NCOLS = (\d+)")
    XY_RE = re.compile(r"F[12][\w]{4,5} = ([-\d.]+)")

    def __init__(self, filename:str) -> None:
        self.filename = filename

    def _get_consts(self, file: FileIO):
        # read the head of the file to find the number of rows, columns and 
        # the limits of the X and Y axes
        header = file.read(600)
        file.seek(0)
        #find rows and columns
        self.NROWS = int(re.findall(self.ROWS_RE, header)[0])
        self.NCOLS = int(re.findall(self.COLS_RE, header)[0])
        #find axis limits
        self.xmax, self.xmin, self.ymax, self.ymin = [float(value) for value in re.findall(self.XY_RE, header)]
        
    def read_file(self, verbose = False):
        if verbose:
            print("Reading Data")
        # reads data and constants from the text file 
        self.data = []
        with open(self.filename, "r") as f:
            self._get_consts(f)
            for line in f:
                if line.startswith("#"):
                    continue
                else:
                    row = [float(line)] + [float(f.readline().strip()) for _ in range(self.NCOLS-1)]
                    self.data.append(row)
        self.data = np.array(self.data)
        print("Shape", self.data.shape)
        print("NCOLS, NROWS", self.NCOLS, self.NROWS)
        
    def process(self):
        x = np.linspace(self.xmax, self.xmin, num=self.NCOLS)
        y = np.linspace(self.ymax, self.ymin, num=self.NROWS)

        arr = self.data
        # expands x and y to the same s
        X, Y = np.meshgrid(x, y)  
        z = arr.flatten()
        # normalise z data
        z = (z/z.max()) * 1
        x = X.flatten()
        y = Y.flatten()
        return x,y,z

def create_mesh(x:np.ndarray, y:np.ndarray, z:np.ndarray) -> mesh.Mesh:
    print("performing Triangulation")
    tri=triang.Triangulation(y,x)
    data = np.zeros(len(tri.triangles), dtype=mesh.Mesh.dtype)
    print("Constructing Mesh")
    NMR_mesh = mesh.Mesh(data, remove_empty_areas=False)
    NMR_mesh.x = x[tri.triangles]
    NMR_mesh.y[:] = y[tri.triangles]
    NMR_mesh.z[:] = z[tri.triangles]
    return NMR_mesh

def main(filename):
    
    spectrum = Bruker2D(filename)
    spectrum.read_file(verbose=True)
    x,y,z = spectrum.process()
    NMR_mesh = create_mesh(x, y, z)
    _, file = os.path.split(filename)
    file, _ = os.path.splitext(file)
    NMR_mesh.save(f'./{file}.stl')

    

if __name__ == "__main__":
    import sys
    try:
        filename = sys.argv[1]
    except:
        filename = "../NMR2.txt"
    main(filename)