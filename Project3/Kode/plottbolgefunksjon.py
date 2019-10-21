from numpy import *
import matplotlib.pyplot as ppl

n=1000
a=2
r3=linspace(-5,5,n)

partikkel1=exp(-a*r3)
partikkel2=exp(a*r3)

ppl.ylim([0, 1])
ppl.plot(r3,partikkel1)
ppl.plot(r3,partikkel2)
ppl.xlabel("Avstand (r)")
ppl.ylabel("Psi^2 (1=100%)")
ppl.show()
print(exp(a*-3))
