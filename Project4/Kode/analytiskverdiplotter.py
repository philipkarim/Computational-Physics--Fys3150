from numpy import *
import matplotlib.pyplot as ppl

n=10000
J=1
k=1
T=linspace(0.8,2.0,n)
b=1/(T*k)


E=-8*J*sinh(8*b*J)/(cosh(8*b*J)+3)
E2=64*J*J*cosh(8*b*J)/(cosh(8*b*J)+3)
M_abs=(2*e**(8*b*J)+4)/(cosh(8*b*J)+3)
M=0
M2=(8*e**(8*b*J)+8)/(cosh(8*b*J)+3)
chi=(M2)/(k*T)
cv=(E2-E**2)/(k*T*T)
    

ppl.grid()
ppl.plot(T,E/4)
ppl.xlabel("Temperature T")
ppl.ylabel(r"$\langle E \rangle$")
ppl.show()

ppl.grid()
ppl.plot(T,M_abs/4)
ppl.xlabel("Temperature T")
ppl.ylabel(r"$\langle |M| \rangle$")
ppl.show()

ppl.grid()
ppl.plot(T,cv/4)
ppl.xlabel("Temperature T")
ppl.ylabel(r"$C_V$")
ppl.show()

ppl.grid()
ppl.plot(T,chi/4)
ppl.xlabel("Temperature T")
ppl.ylabel(r"$\chi$")
ppl.show()

