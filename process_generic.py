import numpy as np
import pandas as pd
import re
from io import FileIO


class Bruker2D:
    ROWS_RE = re.compile(r"NROWS = (\d+)")
    COLS_RE = re.compile(r"NCOLS = (\d+)")
    XY_RE = re.compile(r"F[12][\w]{4,5} = ([-\d.]+)")

    def __init__(self, filename: str) -> None:
        self.filename = filename

    def _get_consts(self, file: FileIO):
        # read the head of the file to find the number of rows, columns and
        # the limits of the X and Y axes
        header = file.read(600)
        file.seek(0)
        # find rows and columns
        self.NROWS = int(re.findall(self.ROWS_RE, header)[0])
        self.NCOLS = int(re.findall(self.COLS_RE, header)[0])
        # find axis limits
        self.xmax, self.xmin, self.ymax, self.ymin = [
            float(value) for value in re.findall(self.XY_RE, header)
        ]

    def read_file(self, verbose=False):
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
                    row = [float(line)] + [
                        float(f.readline().strip()) for _ in range(self.NCOLS - 1)
                    ]
                    self.data.append(row)
        self.data = np.array(self.data)

    def process(self, max_height=1):
        x = np.linspace(self.xmax, self.xmin, num=self.NCOLS)
        y = np.linspace(self.ymax, self.ymin, num=self.NROWS)

        # expands x and y to the same s
        x, y = np.meshgrid(x, y)
        z = self.data.flatten()
        # normalise z data
        z = (z / z.max()) * max_height
        x = x.flatten()
        y = y.flatten()
        return x, y, z


class Generic:
    def __init__(self, filename) -> None:
        self.filename = filename

    def read_file(self, verbose=False):
        if verbose:
            print("Reading file")
        df = pd.read_csv("./DB_Test.txt", delimiter="\t")
        df = df.dropna(thresh=5, axis=0)
        self.data = df

    def process(self, trunc=3, smoothing=False):
        df = self.data
        if smoothing:
            df = df.rolling(window=3, center=True).mean()

        x = df.keys()[0]
        x = df[x].to_numpy()
        y = np.linspace(x.min(), x.max(), len(df.keys()) - 1)
        x, y = np.meshgrid(x, y)

        z = df[df.keys()[1:]].to_numpy()
        x = x.flatten()[::trunc]
        y = y.flatten()[::trunc]
        z = z.flatten()[::trunc]
        z = (z / max(z)) * 1

        return x, y, z
