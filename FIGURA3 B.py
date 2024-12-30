#!/usr/bin/env python
# coding: utf-8

# In[1]:


from numpy import array,arange,cross, pi, sqrt, random, linspace, cos ,sin,linalg, log, dot
import pandas as pd


# In[2]:


#ctes
alpha = 0.005
phi   = 0.3 
h     = 0.1
beta0 = 0.05
w     = 0.08 #No es necesariamente la frecuencia natural del sistema, fue puesta a conveniencia

beta1p   = arange(0.12,0.23+1*10**(-4),1*10**(-4),float)
print(len(beta1p))


# ## Condición inicial magnetización

# In[3]:


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


# ## Condición inicial del deltaM
# 

# In[4]:


PrimeraIteracion = 85000
UltimaIteracion = 90001

#Valor inicial Mi
random.seed(7)
[Mix,Miy,Miz] = 0.5*(-1.0+2.0*random.rand(3))
    
Mi= array([Mix, Miy, Miz])

#Valor inicial M

mx0 = MX #(El valor inicial del tiempo para el cilo)
my0 = MY
mz0 = MZ
M0 = array([mx0,my0,mz0])
    
#Le quito la proyección sobre M y queda en norma 1

deltaM = (Mi-dot(Mi,M0)*M0)/linalg.norm((Mi-dot(Mi,M0)*M0))


# ## deltaM punto

# In[5]:


def dDelta(DeltaM, MX, MY, MZ,t):

    return dot(J(MX, MY, MZ,t), DeltaM)


# ## Crear tabla

# In[6]:


Lyapunov = pd.DataFrame()
Betapoints = []
Lyapoints = []


# ## Integrar

# In[ ]:


#limites integracion 
a = 0.0
b = 1200.0
T = 2*pi/w
dt = T/500
niteraciones = 100 #Cada cuantas iteraciones calculo el exp
PrimeraIteracion = 0
UltimaIteracion = 90000
tp = arange(a,b*T,dt)

#Inicializar variables

#Runge-Kutta

for beta1 in beta1p:
    
    #Inicio el array
    ExpLyapunov = 0

    n=0 #Número de términos que se suman para calcular el exp final
    
    
    #Meto todo el código dentro del cilo
    
    def f(M,t):
    #Ec : ṁ = −m x H - alpha m x (m x H)
    #Con: m y H vectores
    #Con: H = h + beta mz (ez)
    #Con h = h(cos(phi)(ez) + sin(phi)(ex))
    #Con beta= beta0 + beta1 cos(wt)
        mx = M[0]
        my = M[1]
        mz = M[2]
    
        beta=beta0 + round(beta1,4)*cos(w*t)
    
        H = h*array([sin(phi),0,cos(phi)]) + array([0,0, beta*mz])
    
        #fm = ṁ
        fm=cross(-M,H)- alpha*cross(M,(cross(M,H)))
    
        return fm
    
    def J(MX, MY, MZ,t):

        beta=beta0 + beta1*cos(w*t)

        d11 = -MZ*alpha*(MZ*beta + h*cos(phi))
    
        d12 = -MZ*beta - h*cos(phi) + 2*h*MY*alpha*sin(phi)
    
        d13 = -(MY + 2*MX*MZ*alpha)*beta - h*MX*alpha*cos(phi) + 2*h*MZ*alpha*sin(phi)
    
        d21 = MZ*beta + h*cos(phi) - h*MY*alpha*sin(phi)
    
        d22 = -alpha*(MZ**2*beta + h*MZ*cos(phi) + h*MX*sin(phi))
    
        d23 = (MX-2*MY*MZ*alpha)*beta - h*(MY*alpha*cos(phi) + sin(phi))
    
        d31 = alpha*(2*MX*MZ*beta + 2*h*MX*cos(phi) - h*MZ*sin(phi))
    
        d32 = 2*MY*alpha*(MZ*beta + h*cos(phi)) + h*sin(phi)
    
        d33 = (MX**2 + MY**2)*alpha*beta - h*MX*alpha*sin(phi)
    
        return array([[d11,d12,d13],[d21,d22,d23], [d31,d32,d33]])

    it  = arange(PrimeraIteracion,UltimaIteracion,1,int)
    for i in it:
        if i==0 or i>50000:
            
            if (i%niteraciones)==0:
        
                #Le quito la proyección sobre M y se renormaliza
                CoefExp = (1/(dt*niteraciones)) * log(linalg.norm(deltaM-dot(deltaM,m)*m))
        
                ExpLyapunov += CoefExp
                n += 1
        
                deltaM = (deltaM-dot(deltaM,m)*m)/linalg.norm(deltaM-dot(deltaM,m)*m)

        k1 = dt*f(m,tp[i])
        k2 = dt*f(m+0.5*k1,tp[i]+0.5*dt)
        k3 = dt*f(m+0.5*k2,tp[i]+0.5*dt) 
        k4 = dt*f(m+k3,tp[i]+dt)
        m += (k1+2*k2+2*k3+k4)/6
    
        k1delta = dt*dDelta(deltaM,m[0],m[1],m[2],tp[i])
        k2delta = dt*dDelta(deltaM+0.5*k1delta,m[0],m[1],m[2],tp[i]+0.5*dt)
        k3delta = dt*dDelta(deltaM+0.5*k2delta,m[0],m[1],m[2],tp[i]+0.5*dt) 
        k4delta = dt*dDelta(deltaM+k3delta,m[0],m[1],m[2],tp[i]+dt)
        deltaM += (k1delta+2*k2delta+2*k3delta+k4delta)/6      

    #Largest Lyapunov Exponent

    LLE = ExpLyapunov/n
    
    #Armar la tabla 
    
    Betapoints.append(round(beta1,4))
    Lyapoints.append(LLE)

Lyapunov["Exp"] = Lyapoints
Lyapunov["Beta1"] = Betapoints


# In[ ]:


Lyapunov.to_csv("Lyapunov-B.csv")


# ## Graficar

# In[ ]:


from pylab import plot,xlabel,ylabel,show,title, figure,savefig,scatter,legend, xlim, ylim, xticks, yticks

Lyapunov = pd.read_csv("Lyapunov-B.csv")

#Grafico mx
figure(figsize=(13,6))

scatter(Lyapunov["Beta1"],Lyapunov["Exp"], s=5, color= "black")

title("LLE vs ${β}_{1}$",fontsize=20)
xlabel("${β}_{1}$",fontsize=16)
ylabel("LLE",fontsize=16)
    
xlim(0.12,0.22)
ylim(-0.002,0.02)
xticks(fontsize=13)
yticks(fontsize=13)   
    
show()

