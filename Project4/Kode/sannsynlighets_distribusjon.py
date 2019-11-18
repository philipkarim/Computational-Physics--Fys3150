import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt

def readfile(filename):
    
    data = np.loadtxt(filename,unpack=True) 
    
    MC_cycles = data[0]            
    T = data[1]                     
    E_avg = data[2]               
    M_absavg = data[4]              
    C_v = data[3]                   
    X = data[5] 		    
    AC = data[6]		   
    
    return T,MC_cycles,E_avg,M_absavg,C_v,X,AC


L = 20; temp = 1.0; random = 1;

filenameMC = 'ExpectationValues_MC_%d_%.1f_%d.txt' % (L,temp,random)
filenameE = 'Energy_MC_%d_%.1f_%d.txt' % (L,temp,random)

T,MC_cycles,E_avg,M_absavg,C_v,X,AC = readfile(filenameMC)
E = np.loadtxt(filenameE,unpack=True)

results, edges = np.histogram(E, normed=True)
binWidth = edges[1] - edges[0]
plt.bar(edges[:-1], results*binWidth, binWidth,color=(0.0, 0.75, 1.0),edgecolor='black')
plt.xlabel(r'$E$')
plt.ylabel(r'$P(E)$')
plt.savefig(r'Probability_%d_%.2f_%.2f.png'%(L,temp,random))
plt.show()
