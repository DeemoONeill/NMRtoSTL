# -*- coding: utf-8 -*-
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import nmrglue as ng


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
