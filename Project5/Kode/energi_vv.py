import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt;
import numpy as np;

def readFile(filename):
    data = np.loadtxt(filename,unpack=True)

    kinetisk_energi = data[5]
    potensiell_energi = data[6]
    total_energi = data[7]
    angular_moment = data[8]
    return kinetisk_energi,potensiell_energi,total_energi,angular_moment

k,p,t,a = readFile("Earth.txt")
N = 100; FinalTime = 20;
h = FinalTime/N;
tid = np.linspace(0,10,N-1)
print(a)

# plt.plot(tid,k,label = "Kinetisk energi")
# plt.plot(tid,p,label = "Potensiell energi")
# plt.ylabel("Energi [Joule/AU]")
# plt.xlabel("Tid [Ar]")
# plt.legend()
# plt.savefig("energi_e.png")
plt.plot(tid,a)
plt.savefig("angular_moment_vv.png")
plt.show()
