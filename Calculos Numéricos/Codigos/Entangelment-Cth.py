import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

def Msqrt(M): ### Método de Newton-Raphson de la raiz cuadrada.
    x = M;
    e = 1;
    while(e > 10E-8):
        x = 0.5*x+0.5*np.dot(M,LA.inv(x));
        e = LA.norm(np.dot(x,x)-M);
    return x   

#m = 0.1;
n = 10000;
x = np.zeros(n);
for i in range(n):
  x[n-i-1]=np.pi*np.cos((2*i+1)*np.pi/(2*n));

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

# Rmax = 50;
# res1 = [];
# res2 = [];
# for k in range(Rmax):
#     aux1,aux2 = Int(k);
#     res1.append(aux1);
#     res2.append(aux2);

# R = 50;
# S = [];
# for r in range(1,R):
#     X = np.zeros((r,r));
#     P = np.zeros((r,r));
#     for i in range(r):
#         for j in range(i,r):
#             X[i,j] = res1[j-i];
#             X[j,i] = X[i,j];
#             P[i,j] = res2[j-i];
#             P[j,i] = P[i,j];
#     Cv = Msqrt(np.dot(X,P));
#     d1 = LA.eigvals(Cv+0.5*np.eye(len(Cv)));
#     d2 = LA.eigvals(Cv-0.499999*np.eye(len(Cv)));
#     S.append(np.dot(d1,np.log(d1))-np.dot(d2,np.log(d2)));
    
#p,o = np.polyfit(range(2,R),S,1);   
# plt.plot(range(1,len(S)+1),S,'o-')
# plt.show()

#np.savez('EE-cth-0.1',S,R,m,x)
# Masas = np.linspace(0.01,1,10);
# Sm=[];
# for i in range(len(Masas)):
#   m = Masas[i];
#   R = int(10/Masas[i]);
#   res1 = [];
#   res2 = [];
#   for k in range(R):
#     aux1,aux2 = Int(k);
#     res1.append(aux1);
#     res2.append(aux2);
#   X = np.zeros((R,R));
#   P = np.zeros((R,R));
#   for k in range(R):
#     for j in range(k,R):
#       X[k,j] = res1[j-k];
#       X[j,k] = X[k,j];
#       P[k,j] = res2[j-k];
#       P[j,k] = P[k,j];
#   Cv = Msqrt(np.dot(X,P));
#   d1 = LA.eigvals(Cv+0.5*np.eye(len(Cv)));
#   d2 = LA.eigvals(Cv-0.499999*np.eye(len(Cv)));
#   Sm.append(np.dot(d1,np.log(d1))-np.dot(d2,np.log(d2)));
  
# plt.plot(np.log(Masas),Sm,'o-')
# plt.show()
# p,o = np.polyfit(np.log(Masas),Sm,1)
#np.savez('EE-cth-masas',Sm,Masas,x)

m = 10E-7;
Rmax = 100;
res1 = [];
res2 = [];
for k in range(Rmax):
    aux1,aux2 = Int(k);
    res1.append(aux1);
    res2.append(aux2);

# X = np.zeros((Rmax,Rmax));
# P = np.zeros((Rmax,Rmax));

# for i in range(Rmax):
# 	for j in range(i,Rmax):
# 		X[i,j] = res1[j-i];
# 		X[j,i] = X[i,j];
# 		P[i,j] = res2[j-i];
# 		P[j,i] = P[i,j];

plt.plot(res1,'o-')

# S = [];

# for r in range(10,Rmax):
# 	Xv = X[0:r,0:r];
# 	Pv = P[0:r,0:r];
# 	Cv = Msqrt(np.dot(Xv,Pv));
# 	d1 = LA.eigvals(Cv+0.5*np.eye(len(Cv)));
# 	d2 = LA.eigvals(Cv-0.499999*np.eye(len(Cv)));
# 	S.append(np.dot(d1,np.log(d1))-np.dot(d2,np.log(d2)));

# plt.plot(range(10,Rmax),S,'o-')
# plt.show()