import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

m = 1;

def f(x,k):
    res = m*np.cos(k*x)/(4*np.pi*np.sqrt(m**2+np.sin(x)**2));
    return res

def g(x,k):
    res = np.sin(k*x)*np.sin(x)/(4*np.pi*np.sqrt(m**2+np.sin(x)**2));
    return res

x = np.linspace(-np.pi,np.pi,100000);

def Int(k):
    I = 0;
    J = 0;
    for i in range(len(x)-1):
        a = x[i];
        b = x[i+1];
        I = I + ((b-a)/6)*(f(a,k)+4*f((a+b)/2,k)+f(b,k));
        J = J + ((b-a)/6)*(g(a,k)+4*g((a+b)/2,k)+g(b,k));
    return I,J

rini =2;
R = 100;
res1 = [];
res2 = [];
for k in range(R):
    aux1,aux2 = Int(k);
    res1.append(aux1);
    res2.append(aux2);
    
A = np.zeros((R,R));
B = np.zeros((R,R));
for i in range(R):
  for j in range(i,R):
    A[i,j] = res1[j-i];
    A[j,i] = A[i,j];
    B[i,j] = res2[j-i];
    B[j,i] = B[i,j];

#C = np.kron(A,[[0,1],[1,0]])-1j*np.kron(B,[[1,0],[0,-1]]);

S = [];
for r in range(rini,R):
  Av = A[0:r,0:r];
  Bv = B[0:r,0:r];
  Cv = 0.5*np.eye(len(Av)*2)-np.kron(Av,[[0,1],[1,0]])-1j*np.kron(Bv,[[1,0],[0,-1]])
  d1 = LA.eigvals(Cv);
  d2 = LA.eigvals(np.eye(len(Cv))-0.99999999*Cv);
  S.append(np.real(-np.dot(d1,np.log(d1))-np.dot(d2,np.log(d2))))

plt.figure(1)
plt.plot(range(rini,R),S,'--')
plt.show()

C = [];
for i in range(len(S)-1):
  C.append((rini+i+1)*(S[i+1]-S[i]));
  
plt.figure(2)
plt.plot(range(rini+1,R),C,'--')
plt.show()

