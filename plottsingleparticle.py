from numpy import *
import matplotlib.pyplot as ppl

n=1000
r1 = linspace(0,5,n)
r2 = linspace(-5,0,n)

partikkel1= exp(-2*(r1))
partikkel2 = exp(2*(r2))

ppl.plot(r1,partikkel1)
ppl.plot(r2,partikkel2)
ppl.title("Bølgefunksjon")
ppl.xlabel("Avstand (r)")
ppl.ylabel("Sannsynlighet(1=100%)")
ppl.show()
