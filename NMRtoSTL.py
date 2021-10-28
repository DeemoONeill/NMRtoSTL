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
