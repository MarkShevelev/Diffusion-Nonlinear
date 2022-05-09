import numpy as np
import matplotlib.pyplot as pt

I, f = np.loadtxt("./data-two-step.txt", dtype = float, unpack = True)
#x = np.loadtxt("./data-one-step-test.txt", dtype = float, unpack = True)
#y = np.loadtxt("./data-no-step-test.txt", dtype = float, unpack = True)

pt.plot(I, f, color = 'red')
#pt.plot(x, color = 'blue')
#pt.plot(y, color = 'black')
pt.show()