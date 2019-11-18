import matplotlib.pyplot as plt;
import numpy as np;


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

filename20 = 'ExpectationValues_T_%d_%d.txt' % (20,0)
filename40 = 'ExpectationValues_T_%d_%d.txt' % (40,0)
filename60 = 'ExpectationValues_T_%d_%d.txt' % (60,0)
filename80 = 'ExpectationValues_T_%d_%d.txt' % (80,0)
filename100 = 'ExpectationValues_T_%d_%d.txt' % (100,0)

T20,MC_cycles20,E_avg20,M_absavg20,C_v20,X20,AC20 = readfile(filename20);
T40,MC_cycles40,E_avg40,M_absavg40,C_v40,X40,AC40 = readfile(filename40);
T60,MC_cycles60,E_avg60,M_absavg60,C_v60,X60,AC60 = readfile(filename60);
T80,MC_cycles80,E_avg80,M_absavg80,C_v80,X80,AC80 = readfile(filename80);
T100,MC_cycles100,E_avg100,M_absavg100,C_v100,X100,AC100 = readfile(filename100);

#def critical_T(X):
#    for i in range(0,len(X),1):
#        if X[i] == max(X):
#            return(T20[i]);
#            
#TC20 = critical_T(X20)
#TC40 = critical_T(X40)
#TC60 = critical_T(X60)
#TC80 = critical_T(X80)
#TC100 = critical_T(X100)
#
#
#plt.plot(T20,E_avg20, label = "L = 20")
#plt.plot(T40,E_avg40, label = "L = 40")
#plt.plot(T60,E_avg60,label = "L = 60")
#plt.plot(T80,E_avg80,label = "L = 80")
#plt.plot(T100,E_avg100,label = "L = 100")
#plt.ylim(-1.8,-1.1)
#plt.xlabel("T")
#plt.ylabel(r'$\langle E \rangle$')
#plt.legend()
#plt.savefig("T_E_0.05.png")
#plt.show()
#
#plt.plot(T20,M_absavg20, label = "L = 20")
#plt.plot(T40,M_absavg40, label = "L = 40")
#plt.plot(T60,M_absavg60,label = "L = 60")
#plt.plot(T80,M_absavg80,label = "L = 80")
#plt.plot(T100,M_absavg100,label = "L = 100")
#plt.xlabel("T")
#plt.ylabel(r'$\langle |\mathscr{M}| \rangle$')
#plt.legend()
#plt.savefig("T_M_0.05.png")
#plt.show()
#
#
#plt.plot(T20,C_v20, label = "L = 20")
#plt.plot(T40,C_v40, label = "L = 40")
#plt.plot(T60,C_v60,label = "L = 60")
#plt.plot(T80,C_v80,label = "L = 80")
#plt.plot(T100,C_v100,label = "L = 100")
#plt.xlabel("T")
#plt.ylabel(r'$C_v$')
#plt.legend()
#plt.savefig("T_Cv_0.05.png")
#plt.show()
#
#
#plt.plot(TC20,max(X20),'bo')
#plt.plot(TC40,max(X40),'bo')
#plt.plot(TC60,max(X60),'bo')
#plt.plot(TC80,max(X80),'bo')
#plt.plot(TC100,max(X100),'bo')

plt.plot(T20,X20, label = "L = 20")
plt.plot(T40,X40, label = "L = 40")
plt.plot(T60,X60,label = "L = 60")
plt.plot(T80,X80,label = "L = 80")
plt.plot(T100,X100,label = "L = 100")
plt.xlabel("T")
plt.ylabel(r'$\chi$')
plt.legend()
plt.savefig("T_X_0.05.png")
plt.show()

#plt.plot(T20,AC20, label = "L = 20")
#plt.plot(T40,AC40, label = "L = 40")
#plt.plot(T60,AC60,label = "L = 60")
#plt.plot(T80,AC80,label = "L = 80")
#plt.plot(T100,AC100,label = "L = 100")
#plt.xlabel("T")
#plt.ylabel('Akseptert spinn')
#plt.legend()
#plt.savefig("T_AC_0.01.png")
#plt.show()
#
#    