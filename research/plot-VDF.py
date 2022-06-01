import numpy as np
import matplotlib.pyplot as pt

I, f = np.loadtxt("./0.01-0.01-400-0.1-VDF.dat", dtype = float, unpack = True)

pt.plot(I, f, color = 'red')
pt.show()