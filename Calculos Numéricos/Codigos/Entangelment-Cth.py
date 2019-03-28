import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

def Msqrt(M): ### Método de Newton-Raphson de la raiz cuadrada.
    x = M;
    for i in range(50):
        x = 0.5*x+0.5*np.dot(M,LA.inv(x));
    return x   

n = 10000;
m = 0.1;
x = np.arange(-np.pi,np.pi,1/n);

def f(x,k):
    res = np.cos(k*x)/np.sqrt(m**2+2*(1-np.cos(x)));
    return res

def g(x,k):
    res = np.cos(k*x)*np.sqrt(m**2+2*(1-np.cos(x)));
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

res1 = [];
res2 = [];
for k in range(700):
    aux1,aux2 = Int(k);
    res1.append(aux1);
    res2.append(aux2);

R = 300;
S = [];
for r in range(2,R):
    X = np.zeros((r,r));
    P = np.zeros((r,r));
    for i in range(r):
        for j in range(i,r):
            X[i,j] = res1[j-i];
            X[j,i] = X[i,j];
            P[i,j] = res2[j-i];
            P[j,i] = P[i,j];
    Cv = Msqrt(np.dot(X,P));
    d1 = LA.eigvals(Cv+0.5*np.eye(len(Cv)));
    d2 = LA.eigvals(Cv-0.4999999*np.eye(len(Cv)));
    S.append(np.dot(d1,np.log(d1))-np.dot(d2,np.log(d2)));
    
m,b = np.polyfit(range(2,R),S,1);
    
plt.plot(range(2,R),S,'--')
plt.show()

np.savez('EE-cth',res1,res2,S,R,m,x)