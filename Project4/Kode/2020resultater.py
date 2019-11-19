import pylab as ppl

infile10 = open('ExpectationValues_MC_20_1.0_0.txt','r')
infile11 = open('ExpectationValues_MC_20_1.0_1.txt','r')
infile20 = open('ExpectationValues_MC_20_2.4_0.txt','r')
infile21 = open('ExpectationValues_MC_20_2.4_1.txt','r')

#skip the first three lines:
#for i in range(3):
#    infile.readline()

MC_list = []     #Monte Carlo cycles
#Temperatur=1, Ordnet
E_list10 = []      #Expectation energy
Cv_list10 = []     #Heatcapacity
M_list10 = []      #Expectation Magnetisation
X_list10=[]        #Succeptibility
ac_list10 = []     #Accepted spins

#Temperatur=1, Uordnet
E_list11 = []      #Expectation energy
Cv_list11 = []     #Heatcapacity
M_list11 = []      #Expectation Magnetisation
X_list11 = []        #Succeptibility
ac_list11 = []     #Accepted spins

#Temperatur=2.4, Ordnet
E_list20 = []      #Expectation energy
Cv_list20 = []     #Heatcapacity
M_list20 = []      #Expectation Magnetisation
X_list20 = []        #Succeptibility
ac_list20 = []     #Accepted spins

#Temperatur=2.4, Uordnet
E_list21 = []      #Expectation energy
Cv_list21 = []     #Heatcapacity
M_list21 = []      #Expectation Magnetisation
X_list21 = []        #Succeptibility
ac_list21 = []     #Accepted spins

for line in infile10:
    words = line.split()
    MC=float(words[0])
    E = float(words[2])
    Cv = float(words[3])
    M = float(words[4])
    X = float(words[5])
    ac = float(words[6])

    MC_list.append(MC)
    E_list10.append(E)
    Cv_list10.append(Cv)
    M_list10.append(M)
    X_list10.append(X)
    ac_list10.append(ac)

infile10.close()

for line in infile11:
    words = line.split()
    E = float(words[2])
    Cv = float(words[3])
    M = float(words[4])
    X = float(words[5])
    ac = float(words[6])

    E_list11.append(E)
    Cv_list11.append(Cv)
    M_list11.append(M)
    X_list11.append(X)
    ac_list11.append(ac)

infile11.close()

for line in infile20:
    words = line.split()
    E = float(words[2])
    Cv = float(words[3])
    M = float(words[4])
    X = float(words[5])
    ac = float(words[6])

    E_list20.append(E)
    Cv_list20.append(Cv)
    M_list20.append(M)
    X_list20.append(X)
    ac_list20.append(ac)

infile20.close()

for line in infile21:
    words = line.split()
    E = float(words[2])
    Cv = float(words[3])
    M = float(words[4])
    X = float(words[5])
    ac = float(words[6])

    E_list21.append(E)
    Cv_list21.append(Cv)
    M_list21.append(M)
    X_list21.append(X)
    ac_list21.append(ac)

infile21.close()


ppl.xlim([0, 10000])
ppl.ylim([-2,-1])
ppl.plot(MC_list,E_list10, label="Ordnede spinn ved T=1")
ppl.plot(MC_list,E_list11, label= "Uordnede spinn ved T=1")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel(r"$\langle E \rangle $")
ppl.savefig('E1.png', bbox_inches='tight')
ppl.figure()

ppl.xlim([0, 10000])
ppl.ylim([-2,-1])
ppl.plot(MC_list,E_list20, label="Ordnede spinn ved T=2.4")
ppl.plot(MC_list,E_list21, label= "Uordnede spinn ved T=2.4")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel(r"$\langle E \rangle $")
ppl.savefig('E24.png', bbox_inches='tight')
ppl.figure()

ppl.xlim([0, 130000])
ppl.ylim([0,8])
ppl.plot(MC_list,Cv_list10, label="Ordnede spinn ved T=1")
ppl.plot(MC_list,Cv_list11, label= "Uordnede spinn ved T=1")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel(r"$C_v$")
ppl.savefig('Cv1ulik.png', bbox_inches='tight')
ppl.figure()

ppl.xlim([0, 130000])
ppl.ylim([0.6,2])
ppl.plot(MC_list,Cv_list20, label="Ordnede spinn ved T=2.4")
ppl.plot(MC_list,Cv_list21, label= "Uordnede spinn ved T=2.4")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel(r"$C_v$")
ppl.savefig('Cv24ulik.png', bbox_inches='tight')
ppl.figure()

ppl.xlim([0, 130000])
ppl.ylim([0,7.8])
ppl.plot(MC_list,Cv_list10, label="Ordnede spinn ved T=1")
ppl.plot(MC_list,Cv_list11, label= "Uordnede spinn ved T=1")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel(r"$C_v$")
ppl.savefig('Cv1lik.png', bbox_inches='tight')
ppl.figure()

ppl.xlim([0, 130000])
ppl.ylim([0,7.8])
ppl.plot(MC_list,Cv_list20, label="Ordnede spinn ved T=2.4")
ppl.plot(MC_list,Cv_list21, label= "Uordnede spinn ved T=2.4")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel(r"$C_v$")
ppl.savefig('Cv24lik.png', bbox_inches='tight')
ppl.figure()


ppl.xlim([0, 40000])
ppl.ylim([0.1,1.1])
ppl.plot(MC_list,M_list10, label="Ordnede spinn ved T=1")
ppl.plot(MC_list,M_list11, label= "Uordnede spinn ved T=1")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel(r"$\langle |M| \rangle $")
ppl.savefig('M1.png', bbox_inches='tight')
ppl.figure()

ppl.xlim([0, 40000])
ppl.ylim([0.1,1.1])
ppl.plot(MC_list,M_list20, label="Ordnede spinn ved T=2.4")
ppl.plot(MC_list,M_list21, label= "Uordnede spinn ved T=2.4")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel(r"$\langle |M| \rangle $")
ppl.savefig('M24.png', bbox_inches='tight')
ppl.figure()


ppl.xlim([0, 130000])
ppl.ylim([-1,70])
ppl.plot(MC_list,X_list10, label="Ordnede spinn ved T=1")
ppl.plot(MC_list,X_list11, label= "Uordnede spinn ved T=1")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel(r"$\chi$")
ppl.savefig('X1ulik.png', bbox_inches='tight')
ppl.figure()

ppl.xlim([0, 130000])
ppl.ylim([0,12])
ppl.plot(MC_list,X_list20, label="Ordnede spinn ved T=2.4")
ppl.plot(MC_list,X_list21, label= "Uordnede spinn ved T=2.4")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel(r"$\chi$")
ppl.savefig('X24ulik.png', bbox_inches='tight')
ppl.figure()

ppl.xlim([0, 130000])
ppl.ylim([-1,70])
ppl.plot(MC_list,X_list10, label="Ordnede spinn ved T=1")
ppl.plot(MC_list,X_list11, label= "Uordnede spinn ved T=1")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel(r"$\chi$")
ppl.savefig('X1lik.png', bbox_inches='tight')
ppl.figure()

ppl.xlim([0, 130000])
ppl.ylim([-1,70])
ppl.plot(MC_list,X_list20, label="Ordnede spinn ved T=2.4")
ppl.plot(MC_list,X_list21, label= "Uordnede spinn ved T=2.4")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel(r"$\chi$")
ppl.savefig('X24lik.png', bbox_inches='tight')
ppl.figure()


ppl.xlim([0, 12000])
ppl.ylim([0,135])
ppl.plot(MC_list,ac_list10, label="Ordnede spinn ved T=1")
ppl.plot(MC_list,ac_list11, label= "Uordnede spinn ved T=1")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel("Aksepterte spinnkonfigurasjoner")
ppl.savefig('ac1.png', bbox_inches='tight')
ppl.figure()

ppl.xlim([0, 12000])
ppl.ylim([0,135])
ppl.plot(MC_list,ac_list20, label="Ordnede spinn ved T=2.4")
ppl.plot(MC_list,ac_list21, label= "Uordnede spinn ved T=2.4")
ppl.legend()
ppl.xlabel("MC sykluser")
ppl.ylabel("Aksepterte spinnkonfigurasjoner")
ppl.savefig('ac24.png', bbox_inches='tight')
ppl.figure()


