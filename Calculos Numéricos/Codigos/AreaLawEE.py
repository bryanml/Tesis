import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

def Msqrt(M):
    x = M;
    for i in range(50):
        x = 0.5*x+0.5*np.dot(M,LA.inv(x));
    return x   

def f(M,n,l):
    A = np.zeros((M,M));
    D = np.zeros((M,M));
    A[0,0] = (1.5)**2
    D[0,0] = l*(l+1);
    for j in range(1,M):
        A[j,j] = ((j+1.5)**2+(j+.5)**2)/((j+1)**2);
        A[j-1,j] = -((j+0.5)**2)/(j*(j+1));
        A[j,j-1] = A[j-1,j];
        D[j,j] = l*(l+1)/((j+1)**2);
    K = A+D;
    D, U = LA.eig(K);
    Uinv = np.transpose(U);
    P = 0.5000000001*np.dot(U,np.dot(np.sqrt(np.diag(D)),Uinv));
    X = 0.5000000001*np.dot(U,np.dot(np.sqrt(np.diag(1/D)),Uinv));
    Pv = P[0:n,0:n];
    Xv = X[0:n,0:n];
    Cv = Msqrt(np.dot(Xv,Pv));
    d = LA.eigvals(Cv);
    Snl = np.dot(d+0.5*np.ones(n),np.log(d+0.5*np.ones(n)))-np.dot(d-0.5*np.ones(n),np.log(d-0.5*np.ones(n)));
    return Snl


res1 = [];
x1 = [];
M = 100;
for n in range(1,50):
    x1.append(n**2);
    q = 0;
    for l in range(1000):
        q = q + (2*l+1)*f(M,n,l);
    res1.append(q);
    
m,b = np.polyfit(x,res,1);

plt.plot(x,res,'o-')
plt.show()
plt.xlabel(r'r^2',fontsize = 16)
plt.title('R = 100, V esfera de radio r',fontsize = 16)
plt.ylabel(r'$S(\rho_{V})$',fontsize = 16)
plt.text(600,650, r'S(r) $= cr^2$',fontsize = 16)
plt.text(600,600,'c = %1.2f'%(m),fontsize = 16)

m,b = np.polyfit(x,res1,1);