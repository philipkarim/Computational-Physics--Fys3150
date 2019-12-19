import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt;
import numpy as np;

def readFile(filename):
    data = np.loadtxt(filename,unpack=True)

    n = data[0]
    theta = data[1]
    return n,theta

n,theta = readFile("presesjon.txt")
plt.plot(n,theta)
plt.ylabel(r"$\theta_p$")
plt.xlabel("Antall baner")
plt.savefig("theta_p.png")
plt.show()
