import pylab as ppl
import numpy as np

def readFile(readFile): 
    x = 0; y = 0;
    file = open(readFile,"r")
    lines = file.readlines()
    for line in lines:
        value = line.split()
        x = (float(value[0]))
        y = (float(value[1]))

    file.close()
    return x,y;

x = np.zeros(7); 
y = np.zeros(7);

x[0],y[0] = readFile("errorFil1")
x[1],y[1] = readFile("errorFil2")
x[2],y[2] = readFile("errorFil3")
x[3],y[3] = readFile("errorFil4")
x[4],y[4] = readFile("errorFil5")
x[5],y[5] = readFile("errorFil6")
x[6],y[6] = readFile("errorFil7")

ppl.plot(x,y, label = "Gauss special")
ppl.legend()
ppl.xlabel("log10(h)")
ppl.ylabel("log10(Maks relativ error)")
ppl.title("log10 plot av relativ error vs. skrittlengde")
ppl.savefig("Errorplot")