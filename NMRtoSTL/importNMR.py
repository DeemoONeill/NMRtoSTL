# -*- coding: utf-8 -*-
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import nmrglue as ng


class importNMR:
    def __init__(self, filename: str) -> None:
        self.filename = filename

    def read_file(self, verbose=False):
        """
        Import 2D NMR data from file.
        Allowed data types are:
            - Bruker 2D processed data (2rr)
            - NMRpipe 2D processed data (.ft2)

        Parameters
        ----------
        verbose : bool, optional
            If true, additional information will be displayed.
            The default is False.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if verbose:
            print("Reading Data")

        # read data and constants from the text file
        import_functions = [
            self.read_Bruker_proc,
            self.read_pipe
            # self.read_Varian
        ]
        for fun in import_functions:
            try:
                self.data, self.udic = fun(self.filename, verbose)
                break
            except:
                pass
        else:
            raise ValueError("File format not recognised")

        # make sure that the data is 2D
        if self.udic["ndim"] != 2:
            raise ValueError("Data must be 2D!")

        # make sure the the data is frequency domain
        if self.udic[0]["freq"] == False or self.udic[1]["freq"] == False:
            raise ValueError(
                "Data must be frequency domain! "
                "Perform a Fourier transform before processing."
            )

    def make_scales(
        self, verbose=False, f1_min=None, f1_max=None, f2_min=None, f2_max=None
    ):
        """
        Creates scales for the X and Y axes in units of ppm.
        If F1 and F2 limits are set then the spectrum will be
        cropped accordingly.

        Parameters
        ----------
        verbose : bool, optional
            If true, additional information will be displayed.
            The default is False.
        f1_min : float, optional
            Lower limit for F1 (Y) axis in ppm.
            If 'None' then the full spectral range will be used.
        f1_max : TYPE, optional
            Upper limit for F1 (Y) axis in ppm.
            If 'None' then the full spectral range will be used.
        f2_min : TYPE, optional
            Lower limit for F2 (X) axis in ppm.
            If 'None' then the full spectral range will be used.
        f2_max : TYPE, optional
            Upper limit for F2 (X) axis in ppm.
            If 'None' then the full spectral range will be used.

        Returns
        -------
        None.

        """
        # make ppm scales
        uc_f1 = ng.fileiobase.uc_from_udic(self.udic, dim=0)
        self.ppm_f1 = uc_f1.ppm_scale()
        uc_f2 = ng.fileiobase.uc_from_udic(self.udic, dim=1)
        self.ppm_f2 = uc_f2.ppm_scale()

        # set limits
        defaults = uc_f1.ppm_limits() + uc_f2.ppm_limits()
        limits_ppm = [f1_max, f1_min, f2_max, f2_min]
        limits_pts = []
        for n, limit in enumerate(limits_ppm):
            if limit == None:
                limit = defaults[n]
                limits_ppm[n] = limit

        # convert limits from ppm to points
        limits_pts.append(uc_f1.i(limits_ppm[0], "PPM"))
        limits_pts.append(uc_f1.i(limits_ppm[1], "PPM"))
        limits_pts.append(uc_f2.i(limits_ppm[2], "PPM"))
        limits_pts.append(uc_f2.i(limits_ppm[3], "PPM"))

        # apply limits
        self.ppm_f1 = self.ppm_f1[limits_pts[0] : limits_pts[1]]
        self.ppm_f2 = self.ppm_f2[limits_pts[2] : limits_pts[3]]
        self.data = self.data[
            limits_pts[0] : limits_pts[1], limits_pts[2] : limits_pts[3]
        ]

        if verbose:
            print("Plotting data")
            self.plot(limits_ppm)

    def read_Bruker_proc(self, filename, verbose=False):
        """
        Import Bruker processed data.

        Parameters
        ----------
        filename : str
            File path for processed data directory '...\pdata\1'
        verbose : bool, optional
            If true, additional information will be displayed.
            The default is False.

        Returns
        -------
        data : numpy array
            2D NMR data array
        udic : dict
            Dictionary containing experimental parameters

        """
        # read the processed data as numpy array
        dic, data = ng.bruker.read_pdata(filename)

        # create a dictionary containing the experiment parameters
        udic = ng.bruker.guess_udic(dic, data, strip_fake=True)

        # fix the time/frequency axis in the dictionary
        udic[0]["time"], udic[1]["time"] = False, False
        udic[0]["freq"], udic[1]["freq"] = True, True
        if verbose:
            print("Bruker frequency data found")
        return data, udic

    def read_pipe(self, filename, verbose=False):
        """
        Import NMRpipe data.

        Parameters
        ----------
        filename : str
            File path for NMRpipe data file '.ft2'
        verbose : bool, optional
            If true, additional information will be displayed.
            The default is False.

        Returns
        -------
        data : numpy array
            2D NMR data array
        udic : dict
            Dictionary containing experimental parameters

        """
        # read the data as numpy array
        dic, data = ng.pipe.read(filename)

        # create a dictionary containing the experiment parameters
        udic = ng.pipe.guess_udic(dic, data)

        if verbose:
            print("NMRpipe data found")
        return data, udic

    """
    def read_Varian(self, filename, verbose=False):
        # read the data as numpy array
        dic, data = ng.varian.read(filename)

        # create a dictionary containing the experiment parameters
        udic = ng.varian.guess_udic(dic, data)

        if verbose:
            print('Varian/Agilent data found')
        return data, udic
    """

    def plot(self, limits):
        """
        Plot 2D contour spectrum.

        Parameters
        ----------
        limits : list of floats
            X and Y limits for plotting 2D data.
            [f1_max, f1_min, f2_max, f2_min]

        Returns
        -------
        None.

        """
        # calculate contour levels
        cmap = matplotlib.cm.Blues_r  # contour map (colors to use for contours)
        contour_start = 85000  # contour level start value
        contour_num = 25  # number of contour levels
        contour_factor = 1.30  # scaling factor between contour levels
        cl = contour_start * contour_factor ** np.arange(contour_num)

        # plot the data
        plt.figure()
        plt.contour(
            self.data,
            cl,
            cmap=cmap,
            extent=limits,
        )
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()
        plt.show()

    def process(self, max_height: int = 1) -> [np.ndarray, np.ndarray, np.ndarray]:
        # expands x and y to the same s
        x, y = np.meshgrid(self.ppm_f1, self.ppm_f2)
        z = self.data.flatten()[::2]
        # normalise z data
        z = (z / z.max()) * max_height
        x = x.flatten()[::2]
        y = y.flatten()[::2]
        return x, y, z
