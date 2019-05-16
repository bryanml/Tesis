import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

def Msqrt(M): ### Método de Newton-Raphson de la raiz cuadrada.
    x = M;
    e = 1;
    while(e > 1E-7):
        x = 0.5*x+0.5*np.dot(M,LA.inv(x));
        e = LA.norm(np.dot(x,x)-M);
    return x   

x = np.linspace(-np.pi,np.pi,10000);
m = 1E-6;

def f(x,k):
    res = np.cos(k*x)/(4*np.pi*np.sqrt(m**2+2*(1-np.cos(x))));
    return res

def g(x,k):
    res = np.cos(k*x)*np.sqrt(m**2+2*(1-np.cos(x)))/(4*np.pi);
    return res

def Int(k):
    I = 0;
    J = 0;
    for i in range(len(x)-1):
        a = x[i];
        b = x[i+1];
        I = I + ((b-a)/6)*(f(a,k)+4*f((a+b)/2,k)+f(b,k));
        J = J + ((b-a)/6)*(g(a,k)+4*g((a+b)/2,k)+g(b,k));
    return I,J

R = 100;

res1 = [];
res2 = [];
for k in range(R):
    aux1,aux2 = Int(k);
    res1.append(aux1);
    res2.append(aux2);
    
X = np.zeros((R,R));
P = np.zeros((R,R));
for i in range(R):
  for j in range(i,R):
    X[i,j] = res1[j-i];
    X[j,i] = X[i,j];
    P[i,j] = res2[j-i];
    P[j,i] = P[i,j];  

S = [];   
rini = 2;
for r in range(rini,R):
  Xv = X[0:r,0:r];
  Pv = P[0:r,0:r];
  Cv = Msqrt(np.dot(Xv,Pv));
  d1 = LA.eigvals(Cv+0.5*np.eye(len(Cv)));
  d2 = LA.eigvals(Cv-0.49999999*np.eye(len(Cv)));
  S.append(np.real(np.dot(d1,np.log(d1))-np.dot(d2,np.log(d2))));

C = [];
for i in range(len(S)-1):
  C.append((rini+i+1)*(S[i+1]-S[i]));
  
plt.plot(range(rini+1,R),C,'--')
plt.show()
