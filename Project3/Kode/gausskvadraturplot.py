import pylab as ppl

Legrende=[0.626595,0.240305,0.190004,0.0158267]
Laguerre=[0.0327256,0.0155979,0.00873567,0.00531721]
y=[10,15,20,25]

ppl.grid("on")
ppl.plot(y,Legrende, label="Legrende- Relativ feil")
ppl.plot(y,Laguerre, label="Laguerre- Relativ feil")
ppl.xlabel("N")
ppl.ylabel("Relativ feil")
ppl.legend()
ppl.show