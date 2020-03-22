import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import core.PlotMethod as PlotMethod

N=60
M=60

x=np.zeros(shape=(M,N))
y=np.zeros(shape=(M,N)) 

#input dara airfoil
data = pd.read_excel('Airfoil2412 baru.xlsx')
airfoil = data.to_numpy()

length=np.int(len(airfoil)) #panjang BC dan FG
length2=np.int((N-length)/2) #panjang AB CD EF dan HG

theta=np.linspace(3/2*np.pi,1/2*np.pi,length) #definisi sudut

for i in range(length): #mendefinisikan titik dari BC sebagai titik airfoil
    x[0,length2+i]=airfoil[i,0]
    y[0,length2+i]=airfoil[i,1]

for i in range(length): #mendefinisikan titik FG
    x[M-1,length2+i]=10*math.cos(theta[i])
    y[M-1,length2+i]=10*math.sin(theta[i])

for i in range(length2):
    x[0,i]=10*(1-i/length2) #titik AB
    y[0,i]=0
    x[0,length+length2+i]=10*(i+1)/length2 #titik CD
    y[0,length+length2+i]=0
    
    x[M-1,i]=10*(1-i/length2) #titik EF
    y[M-1,i]=-10 
    x[M-1,length+length2+i]=10*(i+1)/length2 #titik GH
    y[M-1,length+length2+i]=10 
 

psi=np.zeros((M-1,N-1))
eta=np.zeros((M-1,N-1))


for i in range(M-1):
    psi[i,0]=(math.exp(i/M)-1)/(math.exp(1)-1)
#for j in range(N-1):
    #eta[0,j]=(math.exp(j/N)-1)/(math.exp(1)-1)

for j in range(N-1):
    x[j,0]=10 #titik AH
    y[j,0]=-10*psi[j,0]
    x[j,N-1]=10 #titik DE
    y[j,N-1]=10*psi[j,0]


#nilai k pada boundary
kB=np.zeros((M-1,1))
kT=np.zeros((M-1,1))
kL=np.zeros((1,N-1))
kR=np.zeros((1,N-1))

kBtotal=0
kTtotal=0
kLtotal=0
kRtotal=0

for i in range(1,M-1): #kBtotal
    R = math.sqrt(((x[i,0]-x[i-1,0])**2)+((y[i,0]-y[i-1,0])**2))
    kBtotal=kBtotal+R
for i in range(1,M-1):
    R = math.sqrt(((x[i,0]-x[i-1,0])**2)+((y[i,0]-y[i-1,0])**2))
    kB[i,0] = kB[i-1,0]+R/kBtotal

for i in range(1,M-1): #kTtotal
    R = math.sqrt(((x[i,M-1]-x[i-1,M-1])**2+(y[i,M-1]-y[i-1,M-1])**2))
    kTtotal=kTtotal+R
for i in range(1,M-1): 
    R = math.sqrt(((x[i,M-1]-x[i-1,M-1])**2+(y[i,M-1]-y[i-1,M-1])**2))
    kT[i,0]=kT[i-1,0]+R/kTtotal
    
for j in range(1,N-1): #kLtotal
    R = math.sqrt(((x[0,j]-x[0,j-1])**2)+((y[0,j]-y[0,j-1])**2))
    kLtotal=kLtotal+R
for j in range(1,N-1): 
    R = math.sqrt(((x[0,j]-x[0,j-1])**2)+((y[0,j]-y[0,j-1])**2))
    kL[0,j] = kL[0,j-1]+R/kLtotal    

for j in range(1,N-1): #kRtotal
    R = math.sqrt(((x[N-1,j]-x[N-1,j-1])**2)+((y[N-1,j]-y[N-1,j-1])**2))
    kRtotal=kRtotal+R
for j in range(1,N-1): 
    R = math.sqrt(((x[i,N-1]-x[i-1,N-1])**2)+((y[i,N-1]-y[i-1,N-1])**2))
    kR[0,j] = kL[0,j-1]+R/kRtotal    

#nilai pada domain
k1=np.zeros((M-1,N-1))
k2=np.zeros((M-1,N-1))

for i in range(M-1):
    for j in range(N-1):
        part1=(1-kL[0,j])*kB[i,0]+kL[0,j]*kT[i,0]
        part2=1-((kT[i,0]-kB[i,0])*(kR[0,j]-kL[0,j]))
        k1[i,j]=part1/part2
         
        part1=(1-kB[i,0])*kL[0,j]+kB[i,0]*kR[0,j]
        part2=1-((kR[0,j]-kL[0,j])*(kT[i,0]-kB[i,0]))
        k2[i,j]=part1/part2

#TFI 
for i in range(1,M-1):
    for j in range(1,N-1):
        U=(1-k1[i,j])*x[0,j]+k1[i,j]*x[M-1,j]
        V=(1-k2[i,j])*x[i,0]+k2[i,j]*x[i,N-1]
        UV=(1-k1[i,j])*(1-k2[i,j])*x[0,0]+k1[i,j]*(1-k2[i,j])*x[M-1,0]+(1-k1[i,j])*k2[i,j]*x[0,N-1]+k1[i,j]*k2[i,j]*(x[M-1,N-1])
        x[i,j]=U+V-UV
        U=(1-k1[i,j])*y[0,j]+k1[i,j]*y[M-1,j]
        V=(1-k2[i,j])*y[i,0]+k2[i,j]*y[i,N-1]
        UV=(1-k1[i,j])*(1-k2[i,j])*y[0,0]+k1[i,j]*(1-k2[i,j])*y[M-1,0]+(1-k1[i,j])*k2[i,j]*y[0,N-1]+k1[i,j]*k2[i,j]*(y[M-1,N-1])
        y[i,j]=U+V-UV
"""
#smoothing grid
alpha=np.zeros((M,N))
beta=np.zeros((M,N))
gamma=np.zeros((M,N))
Rx=np.zeros((M,N))
Ry=np.zeros((M,N))
w=1.5
#errorx=0
#errory=0
lastRvalue=1
iteration=0
error=1
residual=[]


while (error>0.001):
    for i in range(1,M-1):
        for j in range(1,N-1):
            alpha[i,j]=((x[i,j+1]-x[i,j-1])**2+(y[i,j+1]-y[i,j-1])**2)/4
            beta[i,j]=(((x[i+1,j]-x[i-1,j])*(x[i,j+1]-x[i,j-1]))+((y[i+1,j]-y[i-1,j])*(y[i,j+1]-y[i,j-1])))/4
            gamma[i,j]=((x[i+1,j]-x[i-1,j])**2+(y[i+1,j]-y[i-1,j])**2)/4

            A=alpha[i,j]*(x[i+1,j]-2*x[i,j]+x[i-1,j])
            B=(2*beta[i,j]*(x[i+1,j+1]-x[i-1,j+1]-x[i+1,j-1]+x[i-1,j-1]))/4
            C=gamma[i,j]*(x[i,j+1]-(2*x[i,j])+x[i,j-1])
            
            Rx[i,j]=A-B+C
            Ry[i,j]=A-B+C
            
            #update point
            x[i,j]=x[i,j]+(w*(Rx[i,j]/(2*(alpha[i,j]+gamma[i,j]))))
            y[i,j]=y[i,j]+(w*(Ry[i,j]/(2*(alpha[i,j]+gamma[i,j]))))
    
    iteration=iteration+1
    
    currentRvalue=np.sqrt(np.sum(Rx)**2 + np.sum(Ry)**2)
    error=abs(lastRvalue-currentRvalue)
    
    lastRvalue=currentRvalue
"""
plt.scatter(x,y,s=2)
PlotMethod.plotGrid(x, y)
