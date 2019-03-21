import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

def cA(dimA,dimB,dimC):
    aux1 = [];
    aux2 = [];
    for i in range(dimA):
        e = np.transpose(np.eye(1,dimA,i));
        aux = np.kron(e,np.kron(np.eye(dimB),np.eye(dimC)));
        aux1.append(aux);
        aux2.append(np.transpose(aux));
    return aux1 , aux2

def cB(dimA,dimB,dimC):
    aux1 = [];
    aux2 = [];
    for i in range(dimB):
        e = np.transpose(np.eye(1,dimB,i));
        aux = np.kron(np.eye(dimA),np.kron(e,np.eye(dimC)));
        aux1.append(aux);
        aux2.append(np.transpose(aux));
    return aux1 , aux2

def cC(dimA,dimB,dimC):
    aux1 = [];
    aux2 = [];
    for i in range(dimC):
        e = np.transpose(np.eye(1,dimC,i));
        aux = np.kron(np.eye(dimA),np.kron(np.eye(dimB),e));
        aux1.append(aux);
        aux2.append(np.transpose(aux));
    return aux1 , aux2

def TrA(M,dimA,dimB,dimC):
    [EA,TA] = cA(dimA,dimB,dimC);
    res = np.zeros((dimB*dimC,dimB*dimC));
    for i in range (dimA):
        res = res+np.dot(TA[i],np.dot(M,EA[i]));
    return res

def TrB(M,dimA,dimB,dimC):
    [EB,TB] = cB(dimA,dimB,dimC);
    res = np.zeros((dimA*dimC,dimA*dimC));
    for i in range (dimB):
        res = res+np.dot(TB[i],np.dot(M,EB[i]));
    return res

def TrC(M,dimA,dimB,dimC):
    [EC,TC] = cC(dimA,dimB,dimC);
    res = np.zeros((dimA*dimB,dimA*dimB));
    for i in range (dimC):
        res = res+np.dot(TC[i],np.dot(M,EC[i]));
    return res

def DM(dim,aux): ## Genera una matriz densidad
    A = np.random.rand(dim,aux);
    B = np.random.rand(dim,aux);
    C= A+1j*B;
    rho = np.dot(C,np.conjugate(np.transpose(C)))/np.trace(np.dot(C,np.conjugate(np.transpose(C))));
    return rho

def S(M): ## Calcula la entropía de M
    l = LA.eigvals(M);
    S = np.real(-np.dot(l,np.log(l)));
    return S

def f(V):
    M = [];
    for i in range(len(V)):
        M.append(np.dot(V[i],np.conjugate(V)));
    return np.array(M)

##### Strong Subadditivity ######

def SSA(rho,dimA,dimB,dimC): ### Devuelve 1 si vale SSA y 0 si no
    rhoBC = TrA(rho,dimA,dimB,dimC);
    rhoC = TrB(rhoBC,1,dimB,dimC);
    res = S(rhoC)+S(rhoBC)-S(rhoC)-S(rho);
    if res > 0:
        return 1, res;
    else:
        return 0,res;
    
dimA = 15;
dimB = 15;
dimC = 15;
dim = dimA*dimB*dimC;
aux = 4;    
x = [];
y = [];
for i in range(100):
    [n,t] = SSA(DM(dim,aux),dimA,dimB,dimC);
    x.append(n);
    y.append(np.real(t));

plt.plot(x,'o-')
plt.show()

##### Page's Theorem #####
n = 20;
mmax =20; 
y = [];
w = [];
k = 1000;
for m in range(2,mmax):
    q = 0;
    p = 0;
    for i in range(k):
#        V = (np.ones(n*m)-2*np.random.rand(n*m))+1j*(np.ones(n*m)-2*np.random.rand(n*m));
        V = np.random.rand(n*m)+1j*np.random.rand(n*m);
        RPE = V/LA.norm(V); ## random pure state
        PDM = f(RPE); ## RPE pure density matrix
        rhoA = TrB(PDM,n,m,1); ## Matriz densidad reducida
        q = q + S(rhoA)/k;
        p = p + np.real(np.trace(np.dot(rhoA,rhoA)))/k;
    w.append(p);    
    y.append(q);

res = [];
e = [];
for m in range(2,mmax):
    e.append(np.log(m)-m/(2*n));
    q = -(m-1)/(2*n);
    for k in range(n+1,n*m+1):
        q = q+1/k;
    res.append(q);
   
plt.plot(range(2,mmax),y,'o-')
plt.show()
plt.plot(range(2,mmax),res)
plt.show()
plt.title(r'DimA = n = 20 (solo coef. positivos en $\vert \Psi \rangle$)',fontsize = 18)
plt.xlabel('DimB=m',fontsize = 18)
plt.ylabel(r'$S_{n,m}$',fontsize = 18)
plt.legend([r'$\langle S_{A} \rangle$',r'$\sum_{k=n+1}^{nm}\frac{1}{k}-\frac{m-1}{2n}$'],fontsize = 18)


res2=[];
for m in range(2,mmax):
    res2.append((m+n)/(n*m+1))

plt.plot(range(2,mmax),w,'o-')
plt.show()
plt.plot(range(2,mmax),res2)
plt.show()
plt.title(r'DimA = n = 20 (solo coef. positivos en $\vert \Psi \rangle$)',fontsize = 18)
plt.xlabel('DimB=m',fontsize = 18)
plt.ylabel(r'$\langle Tr(\rho_A^2) \rangle$',fontsize = 18)
plt.legend([r'$\langle Tr(\rho_A^2) \rangle$',r'$\frac{m+n}{nm+1}$'],fontsize = 18)

