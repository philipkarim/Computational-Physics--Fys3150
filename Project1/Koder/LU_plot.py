import pylab as ppl
import numpy as np

def exact(x):
    return 1.0 - (1 - np.exp(-10))*x -np.exp(-10*x)
def readFile(readFile): 
    liste = [];
    file = open(readFile,"r")
    lines = file.readlines()
    for line in lines:
        value = line.split()
        liste.append(float(value[0]))
    file.close()
    return liste;

def main():
    v10 = readFile("LU_plot1");
    v100 = readFile("LU_plot2");
    v1000 = readFile("LU_plot3");
    
    x10 = np.linspace(0, 1, len(v10));
    x100 = np.linspace(0, 1, len(v100));
    x1000 = np.linspace(0, 1, len(v1000));
    
    x = np.linspace(0,1,1000);
    ppl.plot(x10,v10, label = "n = 10")
    ppl.plot(x100,v100, label = "n = 100")
    ppl.plot(x1000,v1000, label = "n = 1000")
    ppl.plot(x,exact(x),"r-",label = "Exact")

    ppl.legend()
    ppl.xlabel("x")
    ppl.ylabel("v''(x)")
    ppl.title("LU losning av andregrads differentiallikning med ulike n-verdier ")
    ppl.savefig("LU_plot")


main()
