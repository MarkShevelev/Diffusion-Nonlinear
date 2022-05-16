import numpy as np
import matplotlib.pyplot as pt

#I, f = np.loadtxt("./data.txt", dtype = float, unpack = True)
#x = np.loadtxt("./data-one-step-test.txt", dtype = float, unpack = True)
#y = np.loadtxt("./data-no-step-test.txt", dtype = float, unpack = True)

t, avg = np.loadtxt("avg-vs-time.dat", dtype = float, unpack = True)
pt.plot(t, avg, color = "black")

#pt.plot(I, f, color = 'red')
#pt.plot(x, color = 'blue')
#pt.plot(y, color = 'black')
pt.show()