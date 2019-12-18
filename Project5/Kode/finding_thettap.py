import numpy as np;
import math as m;

def readFile(filename):
    data = np.loadtxt(filename, unpack=True)
    x = data[0]
    y = data[1]
    return x,y

filename = 'presesjon.txt'
x,y = readFile(filename)

X = [abs(number) for number in x]
Y = [abs(number) for number in y]
print(min(X))
print(min(Y))
