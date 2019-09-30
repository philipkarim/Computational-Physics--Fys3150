import pylab as ppl

xJ=[0,7.8,24.36,56.119,117.22]
xA=[0,0,0.001957, 0.0029752,0.005992]
y=[0,50,100,150,200]


ppl.grid("on")
ppl.plot(y,xA, label="Armadillo")
ppl.plot(y,xJ, label="Jacobi")
ppl.title("Algoritmetid med hensyn på N")
ppl.xlabel("N")
ppl.ylabel("tid(s)")
ppl.axis("equal")
ppl.show
