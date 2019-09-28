import pylab as ppl
import numpy as np

def readFile(readFile): 
    liste = [];
    file = open(readFile,"r")
    lines = file.readlines()
    for line in lines:
        value = line.split()
        liste.append(float(value[0]))
    file.close()
    return liste;

def main(fil1,fil2,fil3,fil4,fil5):
    list_10 = readFile(fil1)    
    list_20 = readFile(fil2)
    list_30= readFile(fil3)
    list_40= readFile(fil4)   
    list_50= readFile(fil5)
    
    x = np.linspace(10, 50, 5);

    plot_list_1 = np.zeros(5);
    plot_list_1[0] = list_10[0];
    plot_list_1[1] = list_20[0];
    plot_list_1[2] = list_30[0];
    plot_list_1[3] = list_40[0];
    plot_list_1[4] = list_50[0];
    
    plot_list_2 = np.zeros(5);
    plot_list_2[0] = list_10[1];
    plot_list_2[1] = list_20[1];
    plot_list_2[2] = list_30[1];
    plot_list_2[3] = list_40[1];
    plot_list_2[4] = list_50[1];
    
    plot_list_3 = np.zeros(5);
    plot_list_3[0] = list_10[2];
    plot_list_3[1] = list_20[2];
    plot_list_3[2] = list_30[2];
    plot_list_3[3] = list_40[2];
    plot_list_3[4] = list_50[2];
    
    plot_list_4 = np.zeros(5);
    plot_list_4[0] = list_10[3];
    plot_list_4[1] = list_20[3];
    plot_list_4[2] = list_30[3];
    plot_list_4[3] = list_40[3];
    plot_list_4[4] = list_50[3];
    
    ppl.plot(x,plot_list_1,label = "\u03bb1")
    ppl.plot(x,plot_list_2,label = "\u03bb2")
    ppl.plot(x,plot_list_3,label = "\u03bb3")
    ppl.plot(x,plot_list_4,label = "\u03bb4")

    
    ppl.legend()
    ppl.xlabel("\u03C1_max")
    ppl.ylabel("Relativ feil")
    ppl.title("Relativ feil av \u03bb-ene mot \u03C1_max ved N = 300")
    ppl.savefig("Errorplot_300")


#main("fil_100_10","fil_100_20","fil_100_30","fil_100_40","fil_100_50")
#main("fil_200_10","fil_200_20","fil_200_30","fil_200_40","fil_200_50")
main("fil_300_10","fil_300_20","fil_300_30","fil_300_40","fil_300_50")
