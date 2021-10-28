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


class importNMR:
    def __init__(self, filename: str) -> None:
        self.filename = filename

    def read_file(self, verbose=False):
        if verbose:
            print("Reading Data")
        # read data and constants from the text file
        try:  # read Bruker data
            self.data, self.udic = self.readBruker(self.filename)
        except:
            raise ValueError("File format not recognised")
        # make sure that the data is 2D
        if self.udic["ndim"] != 2:
            raise ValueError("Data must be 2D!")
        # make ppm scales
        uc_f1 = ng.fileiobase.uc_from_udic(self.udic, dim=1)
        self.ppm_f1 = uc_f1.ppm_scale()
        uc_f2 = ng.fileiobase.uc_from_udic(self.udic, dim=0)
        self.ppm_f2 = uc_f2.ppm_scale()
        if verbose:
            print("Plotting data")
            self.plot(uc_f1, uc_f2)

    def readBruker(self, filename):
        # read the processed data as numpy array
        dic, data = ng.bruker.read_pdata(filename)
        # create a dictionary containing the experiment parameters
        udic = ng.bruker.guess_udic(dic, data, strip_fake=True)
        return data, udic

    def plot(self, uc_f1, uc_f2):
        # define plot limits
        ppm_f1_left, ppm_f1_right = uc_f1.ppm_limits()
        ppm_f2_left, ppm_f2_right = uc_f2.ppm_limits()
        # calculate contour levels
        cmap = matplotlib.cm.Blues_r  # contour map (colors to use for contours)
        contour_start = 100000  # contour level start value
        contour_num = 25  # number of contour levels
        contour_factor = 1.20  # scaling factor between contour levels
        cl = contour_start * contour_factor ** np.arange(contour_num)
        # plot the data
        plt.figure()
        plt.contour(
            self.data,
            cl,
            cmap=cmap,
            extent=(ppm_f1_left, ppm_f1_right, ppm_f2_left, ppm_f2_right),
        )
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()
        plt.show()

    def process(self, max_height: int = 1) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        # expands x and y to the same s
        x, y = np.meshgrid(self.ppm_f1, self.ppm_f2)
        z = self.data.flatten()[::2]
        # normalise z data
        z = (z / z.max()) * max_height
        x = x.flatten()[::2]
        y = y.flatten()[::2]
        return x, y, z


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
    spectrum = importNMR(filename)
    spectrum.read_file(verbose=True)
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
        filename = r"Example_data\Bruker_COSY\pdata\1"
    main(filename)
