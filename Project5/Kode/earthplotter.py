import pylab as ppl
from math import *

infile1 = open('euler03n100','r')
infile2 = open('vv03n100','r')
infile3 = open('euler03n100','r')
infile4 = open('vv03n1000','r')
infile5 = open('vvv8pi10000','r')
infile6 = open('vvb103292000','r')
infile7 = open('vvb10332000','r')
infile8 = open('vv03n1000','r')

euler_x = []
euler_y = []
vv_x = []
vv_y = []

euler_x1 = []
euler_y1 = []
vv_x1 = []
vv_y1 = []

euler_x2 = []
euler_y2 = []
vv_x2 = []
vv_y2 = []

euler_x3 = []
euler_y3 = []
vv_x3 = []
vv_y3 = []



for line in infile1:
    words = line.split()
    x=float(words[1])
    y= float(words[2])

    euler_x.append(x)
    euler_y.append(y)
infile1.close()

for line in infile2:
    words = line.split()
    x=float(words[1])
    y= float(words[2])

    vv_x.append(x)
    vv_y.append(y)
infile2.close()

for line in infile3:
    words = line.split()
    x=float(words[1])
    y= float(words[2])

    euler_x1.append(x)
    euler_y1.append(y)
infile3.close()

for line in infile4:
    words = line.split()
    x=float(words[1])
    y= float(words[2])

    vv_x1.append(x)
    vv_y1.append(y)
infile4.close()

for line in infile5:
    words = line.split()
    x=float(words[1])
    y= float(words[2])

    euler_x2.append(x)
    euler_y2.append(y)
infile5.close()

for line in infile6:
    words = line.split()
    x=float(words[1])
    y= float(words[2])

    vv_x2.append(x)
    vv_y2.append(y)
infile6.close()

for line in infile7:
    words = line.split()
    x=float(words[1])
    y= float(words[2])

    euler_x3.append(x)
    euler_y3.append(y)
infile7.close()

for line in infile8:
    words = line.split()
    x=float(words[1])
    y= float(words[2])

    vv_x3.append(x)
    vv_y3.append(y)
infile8.close()
"""
#ppl.xlim([-0.05, 0.05])
#ppl.ylim([-1.05,-0.95])
ppl.axis("equal")
ppl.grid()
ppl.plot(0,0,'ro') 
ppl.plot(euler_x,euler_y, label="N=1000 \n E.F.")
ppl.plot(vv_x,vv_y, label="N=100 \n V.V.")
ppl.plot(euler_x1,euler_y1, label="N=1000 \n E.F.")
ppl.plot(vv_x1,vv_y1, label="N=100 \n V.V.")
ppl.legend(loc="upper right")
ppl.xlabel("x(AU)")
ppl.ylabel("y(AU)")
#ppl.savefig('euler15N1000.png', bbox_inches='tight')
ppl.figure()
"""
ppl.axis("equal")
ppl.grid()
ppl.plot(0,0,'ro') 
ppl.plot(euler_x,euler_y, label="E.F.:N=100")
ppl.plot(vv_x,vv_y, label="V.V: N=100")
ppl.plot(euler_x1,euler_y1, label="E.F:N=1000")
ppl.plot(vv_x1,vv_y1, label="V.V:N=1000")
ppl.legend(loc="lower right")
ppl.xlabel("x(AU)")
ppl.ylabel("y(AU)")
#ppl.savefig('103alle.png', bbox_inches='tight')
ppl.figure()

"""
#ppl.axis("equal")
ppl.xlim([-15, 4])
ppl.ylim([-4,8])
ppl.grid()
ppl.plot(0,0,'ro') 
ppl.plot(euler_x,euler_y, label=r"$v_0$=2$\pi$ $\frac{AU}{yr}$")
ppl.plot(vv_x,vv_y, label=r"$v_0$=7 $\frac{AU}{yr}$")
ppl.plot(euler_x1,euler_y1, label=r"$v_0$=8 $\frac{AU}{yr}$")
ppl.plot(vv_x1,vv_y1, label=r"$v_0$=8.5 $\frac{AU}{yr}$")
ppl.plot(euler_x2,euler_y2, label=r"$v_0$=$\sqrt{8}$ $\pi$ $\frac{AU}{yr}$")
ppl.legend(loc="lower left")
ppl.xlabel("x(AU)")
ppl.ylabel("y(AU)")
#ppl.savefig('unnslippshastighetbeta2.png', bbox_inches='tight')
ppl.figure()

ppl.axis("equal")
#ppl.xlim([-15, 4])
#ppl.ylim([-4,8])
ppl.grid()
ppl.plot(0,0,'ro') 
ppl.plot(euler_x,euler_y, label=r"$\beta$=2")
ppl.plot(vv_x,vv_y, label=r"$\beta$=2.2")
ppl.plot(euler_x1,euler_y1, label=r"$\beta$=2.4")
ppl.plot(vv_x1,vv_y1, label=r"$\beta$=2.6")
ppl.plot(euler_x2,euler_y2, label=r"$\beta$=2.8")
ppl.plot(vv_x2,vv_y2, label=r"$\beta$=2.9")
ppl.plot(euler_x3,euler_y3, label=r"$\beta$=3")
ppl.legend(loc="upper right")
ppl.xlabel("x(AU)")
ppl.ylabel("y(AU)")
#ppl.savefig('betaer100.png', bbox_inches='tight')
ppl.figure()

#Error beregninger
x0=1
y0=0

xe=euler_x[-1]
ye=euler_y[-1]
xv=vv_x[-1]
yv=vv_y[-1]
xe1=euler_x1[-1]
ye1=euler_y1[-1]
xv1=vv_x1[-1]
yv1=vv_y1[-1]

xed=abs(x0-xe)
yed=abs(y0-ye)
xed1=abs(x0-xe1)
yed1=abs(y0-ye1)

xvd=abs(x0-xv)
yvd=abs(y0-yv)
xvd1=abs(x0-xv1)
yvd1=abs(y0-yv1)

error_ef=sqrt(xed**2+yed**2)
error_vv=sqrt(xvd**2+yvd**2)

error_ef1=sqrt(xed1**2+yed1**2)
error_vv1=sqrt(xvd1**2+yvd1**2)
error_ef2=sqrt(xed1**2+yed1**2)
error_vv2=sqrt(xvd1**2+yvd1**2)
print("Feil for EF, N=1000", error_ef)
print("Feil for VV, N=1000", error_vv)
print("Feil for EF, N=10000", error_ef1)
print("Feil for VV, N=10000", error_vv1)
print("Feil for EF, N=100000", error_ef2)
print("Feil for VV, N=100000", error_vv2)
"""

