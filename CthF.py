import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
import time

print('m:')

m = input()

print('R:')

R = input()

m = float(m);
R = int(R); 
t0 = time.clock();
def Msqrt(M): ### Método de Newton-Raphson de la raiz cuadrada.
    x = M;
    e = 1;
    while(e > 1E-7):
        x = 0.5*x+0.5*np.dot(M,LA.inv(x));
        e = LA.norm(np.dot(x,x)-M);
    return x   

x = np.linspace(-np.pi,np.pi,100000); ## discretización del intervalo [-pi,pi]
####m = 1E-4;### Masa del bosón

def f(x,k): #### función a integrar para calcular X 
    res = m*np.cos(k*x)/(4*np.pi*np.sqrt(m**2+np.sin(x)**2)); 
    return res

def g(x,k): #### función a integrar para calcular P 
    res = np.sin(x)*np.sin(k*x)/(4*np.pi*np.sqrt(m**2+np.sin(x)**2));
    return res

def Int(k): #### Integral para k = i-j.
    I = 0;
    J = 0;
    for i in range(len(x)-1):
        a = x[i];
        b = x[i+1];
        I = I + ((b-a)/6)*(f(a,k)+4*f((a+b)/2,k)+f(b,k));
        J = J + ((b-a)/6)*(g(a,k)+4*g((a+b)/2,k)+g(b,k));
    return I,J

###R = 50; #### Radio máximo. El sistema es [0,1,2,...,R-1].

res1 = [];
res2 = [];
for k in range(R): #### Integral de todos loss f y g desde 0 hasta R.
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

S = []; ##### Vector de entropias S(r) para un subsistema V = [0,1,...,r]
rini = 2; #### r inicial
for r in range(rini,R):
  Av = A[0:r,0:r];
  Bv = B[0:r,0:r];
  Cv = 0.5*np.eye(2*len(Av))-np.kron([[0,1],[1,0]],Av)-1j*np.kron([[1,0],[0,-1]],Bv);
  d1 = LA.eigvals(Cv); #### autovalores de Cv + 1/2
  d2 = LA.eigvals(np.eye(len(Cv))-0.999999*Cv); #### autovalores de Cv - 1/2
  S.append(-np.real(np.dot(d1,np.log(d1))+np.dot(d2,np.log(d2)))); #### S(r) = Tr[(C+1/2)log(C+1/2)-(C-1/2)log(C-1/2)]

C = []; #### función c = rdS/dr
for i in range(len(S)-1):
  C.append((rini+i+1)*(S[i+1]-S[i]));

t = time.clock()-t0;

np.save('data', [m,R,S,C])

plt.plot(range(rini+1,R),C,'--')
plt.title('t = %1.4f'%(t),fontsize = 16)
plt.show()
