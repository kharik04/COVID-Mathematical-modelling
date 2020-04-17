import numpy as np
import matplotlib.pyplot as plt


class Stack:
    def __init__(self, size):
        self.stackArr = [0]*size
    def add(self, value):
        self.stackArr.pop(0)
        self.stackArr.append(value)
    def get(self):
        return self.stackArr[0]
#for ii,r in enumerate(rate):
r=6000000
population = 5e6
N=population

I=1
S=population-I
Q=0
R=0
Iarr=[I]
Sarr=[S]
Qarr=[Q]
Rarr=[R]
    #r=1000*24*60#interactions/day
p=0.5
tf=1
dt=1/(24*6000)
recovery = 1#delay between desease and recovery in days
incubation = 1
S_i = Stack(int(incubation/dt))
I_i = Stack(int(incubation/dt))
S_t = Stack(int((incubation+recovery)/dt))
I_t = Stack(int((incubation+recovery)/dt))

for i in range(20):
    t=0
    while t<tf:
        S=S+(-2*r*p*I*S/((N-Q)**2))*dt
        I=I+(2*r*p*(I*S-I_i.get()*S_i.get())/((N-Q)**2))*dt
        Q=Q+(2*r*p*(I_i.get()*S_i.get()-I_t.get()*S_t.get())/((N-Q)**2))*dt
        R=R+(2*r*p*(I_t.get()*S_t.get())/((N-Q)**2))*dt



        S_i.add(S)
        I_i.add(I)
        S_t.add(S)
        I_t.add(I)
        t+=dt
    Iarr.append(I)
    Sarr.append(S)
    Qarr.append(Q)
    Rarr.append(R)
    print(i)
plt.plot(Sarr)
plt.plot(Iarr)
plt.plot(Qarr)
plt.plot([0])
plt.xlabel('t')
#plt.ylabel('I(t)')
plt.title('SIR model, p=0.5, r=600000')
plt.show()
