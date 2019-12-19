import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt;
import numpy as np;

def readFile(filename):
    data = np.loadtxt(filename,unpack=True)

    kinetisk_energi = data[1]
    potensiell_energi = data[2]
    total_energi = data[3]
    angular_moment = data[4]
    return kinetisk_energi,potensiell_energi,total_energi,angular_moment

k,p,t,a = readFile("Energi.txt")
N = 10000; FinalTime = 20;
h = FinalTime/N;

k_list = []
p_list = []

for i in range(0,len(k),1):
    if k[i] > 0.003:
        k_list.append(k[i])
    if p[i] < -0.007:
        p_list.append(p[i])


total = []
for i in range(0,len(k_list),1):
    total.append(p_list[i] + k_list[i])
a_list = []
for i in range(0,len(a),1):
    if a[i] > 30:
        a_list.append(a[i])
tid = np.linspace(0,200,len(a_list))

# plt.plot(tid,k_list,label = "Kinetisk energi")
# plt.plot(tid,p_list,label = "Potensiell energi")
# plt.plot(tid,total,label = "Total energi" )
plt.plot(tid,a_list)
plt.ylabel("Angular moment []")
plt.xlabel("Tid [Ar]")
plt.savefig("angular_moment_helesystemet.png")
plt.show()
