
from  matplotlib import pyplot as plt
import numpy as np
#Function for initialization of parameters
def initialize():
    RMin = 0.0
    RMax = 10.0
    Dim = 400
    return RMin, RMax, Dim
# Here we set up the harmonic oscillator potential
def potential(r,h,w):
    r = r + (i + 1) * h;
    return (r* r * w * w) + (1 / r)
def plotting(w):
    #Get the boundary, orbital momentum and number of integration points
    RMin, RMax, Dim = initialize()
    
    #Initialize constants
    Step    = RMax/(Dim+1)
    DiagConst = 2.0 / (Step*Step)
    NondiagConst =  -1.0 / (Step*Step)
    
    #Calculate array of potential values
    v = np.zeros(Dim)
    r = np.linspace(RMin,RMax,Dim)
    for i in range(Dim):
        r[i] = RMin + (i+1) * Step;
        v[i] = potential(r[i],Step,w)
    
    #Setting up a tridiagonal matrix and finding eigenvectors and eigenvalues
    Hamiltonian = np.zeros((Dim,Dim))
    Hamiltonian[0,0] = DiagConst + v[0];
    Hamiltonian[0,1] = NondiagConst;
    for i in range(1,Dim-1):
        Hamiltonian[i,i-1]  = NondiagConst;
        Hamiltonian[i,i]    = DiagConst + v[i];
        Hamiltonian[i,i+1]  = NondiagConst;
    Hamiltonian[Dim-1,Dim-2] = NondiagConst;
    Hamiltonian[Dim-1,Dim-1] = DiagConst + v[Dim-1];
    # diagonalize and obtain eigenvalues, not necessarily sorted
    EigValues, EigVectors = np.linalg.eig(Hamiltonian)
    # sort eigenvectors and eigenvalues
    permute = EigValues.argsort()
    EigValues = EigValues[permute]
    EigVectors = EigVectors[:,permute]

    FirstEigvector = EigVectors[:,0]
    return r,FirstEigvector
r1,liste1 = plotting(0.01)
r2,liste2 = plotting(0.5)
r3,liste3 = plotting(1)
r4,liste4 = plotting(5)

plt.plot(r1, liste1**2,label = r"$\omega_r= %.2f$" %0.01)
plt.plot(r2, liste2**2,label = r"$\omega_r= %.1f$" %0.5)
plt.plot(r3, liste3**2,label = r"$\omega_r= %.1f$" %1)
plt.plot(r4, liste4**2,label = r"$\omega_r= %.1f$" %5)

plt.axis([0,9,0.0, 0.05])
plt.xlabel('\u03C1')
plt.ylabel('$u(\u03C1)^2$')
plt.title('Sannsynlighetestettheten for ingen interaksjon')
plt.legend()
plt.savefig('to_elektroner')
plt.show()
