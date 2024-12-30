from math import cos, sin
from numpy import array,arange,cross, pi, sqrt, random, linspace
import pandas as pd

#ctes
alpha = 0.005
phi   = 0.3 
h     = 0.1
beta0 = 0.05
w     = 0.08 #No es necesariamente la frecuencia natural del sistema, fue puesta a conveniencia

#Para hacer el ciclo:
beta1p   = arange(0.,0.123+1*10**(-4),1*10**(-4),float)

x = pd.DataFrame()
y = pd.DataFrame()
z = pd.DataFrame()

#Para números aleatorios

random.seed(7)
random.rand(2)
random.rand(2)
random.rand(2)
random.rand(2)

[MX,MY] = 0.5*(-1.0+2.0*random.rand(2))
MZ = sqrt(abs(1-(MX**2+MY**2)))

m = array([MX,MY,MZ])
print(m)

#Ahora el ciclo
for beta1 in reversed(beta1p):
    #Meto todo el código dentro del cilo
    
    #Ec : ṁ = −m x H - alpha m x (m x H)
    #Con: m y H vectores
    #Con: H = h + beta mz (ez)
    #Con h = h(cos(phi)(ez) + sin(phi)(ex))
    #Con beta= beta0 + beta1 cos(wt)



    def f(M,t):
    
        mx = M[0]
        my = M[1]
        mz = M[2]
    
        beta=beta0 + beta1*cos(w*t)
    
        H = h*array([sin(phi),0,cos(phi)]) + array([0,0, beta*mz])
    
        #fm = ṁ
        fm=cross(-M,H)- alpha*cross(M,(cross(M,H)))
    
        return fm
    
   #limites integracion 
    a = 0.0
    b = 1200.0
    T = 2*pi/w
    dt = T/500
    N = int((b*T-a)/dt)
    S=int(round(T/dt,0))
     
    #Runge-Kutta
    tp = arange(a,b*T,dt)
    tpoints = []
    xpoints = []
    ypoints = []
    zpoints = []
        

    for t in tp:         
        if t>80000 or t==0:
        
            tpoints.append(t)
            
            xpoints.append(m[0])
            ypoints.append(m[1])
            zpoints.append(m[2])
            

        k1 = dt*f(m,t)
        k2 = dt*f(m+0.5*k1,t+0.5*dt)
        k3 = dt*f(m+0.5*k2,t+0.5*dt) 
        k4 = dt*f(m+k3,t+dt)
        m += (k1+2*k2+2*k3+k4)/6
        

    #Armar la tabla 
    
    #Tiempo
    x["Tiempo"]=tpoints
    y["Tiempo"]=tpoints
    z["Tiempo"]=tpoints
    
    #m
    
    if (round(beta1*10**(6),0)%500)==0:
    
        x[beta1]=xpoints
        y[beta1]=ypoints
        z[beta1]=zpoints


# x.to_csv(f"CI(5){[MX,MY,MZ]}-beta1-mx-0-0.123.csv")
# y.to_csv(f"CI(5){[MX,MY,MZ]}-beta1-my-0-0.123.csv")
# z.to_csv(f"CI(5){[MX,MY,MZ]}-beta1-mz-0-0.123.csv")

#Leer las tabla
x = pd.read_csv(f"CI(5){[MX,MY,MZ]}-beta1-mx-0-0.123.csv")
y = pd.read_csv(f"CI(5){[MX,MY,MZ]}-beta1-my-0-0.123.csv")
z = pd.read_csv(f"CI(5){[MX,MY,MZ]}-beta1-mz-0-0.123.csv")

#Borrar lo que no me sirve e incomoda
x1 = x.drop(['Tiempo'],axis=1)
x2 = x1.drop(['Unnamed: 0'],axis=1)

#Voltear la tabla
mx = x2.T

#Hacer lo mismo con my y mz
y1 = y.drop(['Tiempo'],axis=1)
y2 = y1.drop(['Unnamed: 0'],axis=1)
my = y2.T

z1 = z.drop(['Tiempo'],axis=1)
z2 = z1.drop(['Unnamed: 0'],axis=1)
mz = z2.T

betapoints = []
for beta1 in beta1p:
    if (round(beta1*10**(6),0)%500)==0:
        betapoints.append(beta1)

from pylab import plot,xlabel,ylabel,show,title, figure,savefig,scatter,legend, xlim, ylim, xticks, yticks

#Para graficar
it     = arange(70350,90705,500,int)
#Grafico mx
figure(figsize=(13,3))
for i in it:
    scatter(list(reversed(betapoints)),mx[i], s=5)
    title("${m}_{x}$ vs ${β}_{1}$",fontsize=20)
    xlabel("${β}_{1}$",fontsize=16)
    ylabel("${m}_{x}$",fontsize=16)
    

    xticks(fontsize=13)
    yticks(fontsize=13)   
    
show()
