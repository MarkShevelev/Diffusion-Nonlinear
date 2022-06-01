import numpy as np
import matplotlib.pyplot as pt

t, R = np.loadtxt("./0.01-0.01-400-0.01-R.dat", dtype = float, unpack = True)

pt.plot(t, R, color = 'red')
pt.show()