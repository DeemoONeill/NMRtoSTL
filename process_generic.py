# %%

import numpy as np
import pandas as pd

class Generic:
    
    def __init__(self, filename) -> None:
        self.filename = filename

    def read_file(self, verbose = False):
        if verbose:
            print("Reading file")
        df = pd.read_csv("./DB_Test.txt", delimiter="\t")
        df = df.dropna(thresh=5, axis=0)
        self.data = df

    def process(self, trunc = 3, smoothing = False):
        df = self.data
        if smoothing:
            df = df.rolling(window=3, center=True).mean()
        
        x = df.keys()[0]
        x = df[x].to_numpy()
        y = np.linspace(x.min(), x.max(),len(df.keys())-1)
        x,y = np.meshgrid(x,y)

        z = df[df.keys()[1:]].to_numpy()
        x = x.flatten()[::trunc]
        y = y.flatten()[::trunc]
        z = z.flatten()[::trunc]
        z = (z/max(z))*1

        return x,y,z

