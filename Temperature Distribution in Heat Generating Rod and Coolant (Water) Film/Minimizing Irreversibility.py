from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np 
import math
def my_function(Z):
    X=Z[0]
    Y=Z[1]
    C=Z[2]
    D=Z[3]
    F=np.empty(4)
    F[0]=T-X-Y #X=c1,Y=c2
    F[1]=T-C-D #C=c1',D=c2'
    F[2]=(E*c*(X+(Y*math.exp(-H*L))+(I*L)-T))-(h1*d*((2*k*h2*(((X-C)*L)-(math.exp(-H*L)*(Y-D)/H)))-(q*t*t*h2*L)-(2*k*q*t*L))/((2*k*(h1-h2))+(2*h1*h2*t)))+(G*c*(C+(D*math.exp(-H*L))+(I*L)-T))-(h2*d*((2*k*h1*(((X-C)*L)-((math.exp(-H*L)-1)*(Y-D)/H)))-(q*t*t*h1*L)-(2*k*q*t*L))/((2*k*(h1-h2))+(2*h1*h2*t)))
    # F[3]=(G*c*(C+(D*math.exp(-H*L))+(I*L)-T))-(h2*d*((2*k*h1*(((X-C)*L)-((math.exp(-H*L)-1)*(Y-D)/H)))-(q*t*t*h1*L)-(2*k*q*t*L))/((2*k*(h1-h2))+(2*h1*h2*t)))
    F[3]=(E*c*(X+(Y*math.exp(-H*L))+(I*L)-T))+(G*c*(C+(D*math.exp(-H*L))+(I*L)-T))-q
    return F
d1=40 # float(input('d in mm ='))
t1=10
# float(input('t in mm= ')) 
# v=0.044*10**6
b1=550
# float(input('b in mm= '))
# d1=(v/(b1*2))
k=23
# float(input('thermal conductivity= '))
T=298
# float(input('ti= '))
m=10
# float(input('m dot= '))
L=2
# float(input('length= '))
q=400*10**3
# float(input('q= '))
c=4184
# float(input('c='))
t=0.001*t1
d=0.001*d1
b=0.001*b1
Tmax=0.0
array=[]
a1=b/2-0.001
while(a1>=d+(t/2)):
    array.append(a1)
    a1=a1-0.001
a1=(b/2)+0.001
while(a1<=b-d-(t/2)):
    array.append(a1)
    a1=a1+0.001
array.sort()
Y=[0]
Tsmaxi=0.0
Tsmaxo=0.0
y=-t/2
while(y<=(t/2)):
    Y.append(y)
    y=y+0.001
To1L=[]
To2L=[]
TSi=[]
TSo=[]
# print(array)
n=int(len(array)/2)
mid1=int(n-1)
mid2=n-1
Tmi=2500
Tmo=2500
sgenmax=0.0
SGEN=[]
# v=d1+t1/2-1 #start value
# while (v<=(b1-d1-t1/2)-1): #stop value
#     v=v+1
#     array.append(v) 
#array=np.linspace(t1/2+1,((b1-t1/2)-1),((b1-t1-2)+1))
# arr=[i*0.001 for i in array]
for a in array:
    alp1=d/(a-t/2)
    alp2=d/(b-a-t/2)
    nu1=((math.log(d/(a-t/2),math.exp(1))+2.6)/0.97)
    if(nu1>5.384):
        nu1=5.384
    nu2=((math.log(d/(b-a-t/2),math.exp(1))+2.6)/0.97)
    if(nu2>5.384):
        nu2=5.384
    di1=4*d*(a-t/2)/(2*(d+(a-t/2)))
    di2=4*d*(b-(a-t/2))/(2*(d+(b-a-t/2)))
    f1=24*(1-1.3553*alp1+1.9467*alp1**2-1.7012*alp1**3+0.9564*alp1**4-0.2537*alp1**5)#friction factors
    f2=24*(1-1.3553*alp2+1.9467*alp2**2-1.7012*alp2**3+0.9564*alp2**4-0.2537*alp2**5)
    h1=(nu1*k/di1)
    h2=(nu2*k/di2)
    M=-(((h1*d*q*t*t*h2)+(2*k*q*t*h1*d))/((2*k*(h1-h2))+(2*h1*h2*t)))#k1
    N=(((h1*d*q*t*t*h2)-(2*k*q*t*h2*d))/((2*k*(h1-h2))+(2*h1*h2*t)))#k2
    J=(((2*k*h1*h2*d)/((2*k*(h1-h2))+(2*h1*h2*t))))#alpha
    A=(f2/f1)*((a-(0.5*t))/(d+a-(0.5*t)))*((d+b-a-(0.5*t))/(b-a-(0.5*t)))
    E=((m*(math.sqrt((f2*(a-(0.5*t))*(d+b-a-(0.5*t)))/(f1*(d+a-(0.5*t))*(b-a-(0.5*t)))))*(a-(0.5*t))/(((math.sqrt((f2*(a-(0.5*t))*(d+b-a-(0.5*t)))/(f1*(d+a-(0.5*t))*(b-a-(0.5*t)))))*(a-(0.5*t)))+(b-a-(0.5*t))))) #m1
    G=m*(b-a-(0.5*t))/((math.sqrt(A)*(a-(0.5*t)))+(b-a-(0.5*t)))#m2 
    H=(((E-G)*J)/(E*G*c)) #a
    I=((M-N)/((c*E)-(c*G))) #b
    ZGuess=np.array([100,1,100,1]) 
    Z = fsolve(my_function,ZGuess)
    To1=Z[0]+Z[1]*math.exp(-H*L)+I*L
    To2=Z[2]+Z[3]*math.exp(-H*L)+I*L
    # if a>(b/2):
    #     To1=To2L[mid1]
    #     To2=To1L[mid1]
    #    mid1=mid1-1
    P=((q*t/2)*(h1+h2))/(k*(h1-h2)+(t*h1*h2)) #c1
    Q=(((h1*t/2)-k)*((q*t/2)+h2*((q*t**2/8*k)+T))-((h2*t/2)+k)*((q*t/2)-h1*((q*t**2/8*k)+T)))/(k*(h1-h2))+h1*h2*t #c2
    R=(((q*t/2)*(h1+h2))+h1*h2*(To2-To1))/(k*(h1-h2)+(t*h1*h2))
    S=(((h1*t/2)-k)*((q*t/2)+h2*((q*t**2/8*k)+To2))-((h2*t/2)+k)*((q*t/2)-h1*((q*t**2/8*k)+To1)))/(k*(h1-h2))+h1*h2*t
    TS=[]
    for y in Y:
        Tsi=(-q/k)*0.5*y**2+P*y+Q
        Tso=(-q/k)*0.5*y**2+R*y+S
        if (y==-(t/2) or y==(t/2)):
              TS.append(Tsi)
              TS.append(Tso)      
        if (Tsi>Tsmaxi):
            Tsmaxi=Tsi
            maxposi=y
        if (Tso>Tsmaxo):
            Tsmaxo=Tso
            maxposo=y
    if a>(b/2):
        Tsmaxi=TSi[mid2]
        Tsmaxo=TSo[mid2]
        mid2=mid2-1
    if a==(b/2)+0.001:
        # print(Z,H,I)
        TM=[Tsmaxi,Tsmaxo]
    if a<(b/2):
        if Tsmaxi<Tmi:
            Tmi=Tsmaxi
            ai=a
            Zi=Z
            Hi=H
            Ii=I
        if Tsmaxo<Tmo:
            Tmo=Tsmaxo
            ao=a
            Zo=Z
            Ho=H
            Io=I
            TB=[To1,To2]
    TSi.append(Tsmaxi)
    TSo.append(Tsmaxo)
    To1L.append(To1)
    To2L.append(To2)
    x=0.5*10**(-3)
    sgen1=0.0
    sgen2=0.0
    while(x<L):
        Tb1=Z[0]+Z[1]*math.exp(-H*x)+I*x
        Tb2=Z[2]+Z[3]*math.exp(-H*x)+I*x
        c1=(((q*t/2)*(h1-h2))+h1*h2*(Tb2-Tb1))/(k*(h1+h2)+(t*h1*h2))
        c2=((h1*t/2)+k)*((q*t/2)+h2*((q*t**2/8*k)+Tb2))+((h2*t/2)+k)*((q*t/2)+h1*((q*t**2/8*k)+T))/((k*(h1-h2))+h1*h2*t)
        Ts1=(-q/k)*0.5*(-t/2)**2+c1*(-t/2)+c2
        Ts2=(-q/k)*0.5*(t/2)**2+c1*(t/2)+c2
        sgentemp1=h1*(h2*q*t**2+2*k*q*t-2*k*h2*(Tb1-Tb2)/(2*k*(h1+h2)+2*h1*h2*t))*d*0.001*((1/Tb1)-(1/Ts1))-E*f1*0.001*(E/(1000*(a-0.5*t)*d)**2/2*9.81*di1*1000)
        sgentemp2=h2*(h1*q*t**2+2*k*q*t-2*k*h1*(Tb2-Tb1)/(2*k*(h1+h2)+2*h1*h2*t))*d*0.001*((1/Tb2)-(1/Ts2))-G*f2*0.001*(E/(1000*(b-a-0.5*t)*d)**2/2*9.81*di2*1000)
        sgen1+=sgentemp1
        sgen2+=sgentemp2
        x+=0.001
    sgen=sgen1+sgen2
    SGEN.append(sgen)
    if (sgen>sgenmax):
        sgenmax=sgen
        maxp=a
    print(To1,To2,end="\t")
    print(a,end="\t")
    # print(E,G,end="\t")
    print(Tsmaxi,maxposi,end="\t")
    print(Tsmaxo,maxposo,end="\t")
    print(TS)
print(ai,ao,Tmi,Tmo)
print(TM,(b1/d1))
print(TB)
print(sgenmax,maxp)
print(SGEN)
plt.plot(array,SGEN)
#plt.plot(array,SGEN,'o')
plt.show()
# print("integration const at aopt",Zo,"a=",Ho,"b=",Io)