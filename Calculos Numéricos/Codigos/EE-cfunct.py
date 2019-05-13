import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

x = np.arange(-np.pi,np.pi,0.0001);

def fX(x,k):
    res = np.cos(k*x)/(4*np.pi*np.sqrt(m**2+2*(1-np.cos(x))));
    return res

def fP(x,k):
    res = np.cos(k*x)*np.sqrt(m**2+2*(1-np.cos(x)))/(4*np.pi);
    return res

def IX(x,k):
    res = 0;
    for i in range(len(x)-1):
        res = res + (fX(x[i],k)+4*fX((x[i]+x[i+1])/2,k)+fX(x[i+1],k))*(x[i+1]-x[i])/6;
    return res

def IP(x,k):
    res = 0;
    for i in range(len(x)-1):
        res = res + (fP(x[i],k)+4*fP((x[i]+x[i+1])/2,k)+fP(x[i+1],k))*(x[i+1]-x[i])/6;
    return res

def Msqrt(M):
    x = M;
    for i in range(50):
        x = 0.5*x+0.5*np.dot(M,LA.inv(x));
    return x   

##### X y P en función de la masa #####
#masas =[0.001,0.01,0.1,1,10];
#leg = [];
#X = [];
#P = [];
#for l in range(len(masas)):
#    aux1 = [];
#    aux2 = [];
#    m = masas[l];
#    leg.append(str(m));
#    for k in range(40):
#        aux1.append(IX(x,k));
#        aux2.append(IP(x,k));
#    X.append(aux1);
#    P.append(aux2);
#    
#plt.figure(1)
#
#for l in range(len(masas)):
#    plt.plot(X[l],'o-')
#    plt.show()
#    
#plt.legend(leg)  
#    
#plt.figure(2)
#
#for l in range(len(masas)):
#    plt.plot(P[l],'o-')
#    plt.show()
#    
#plt.legend(leg) 

#### a) #######
#m = 0.01;
#X = [];
#P = [];
#for k in range(500):
#        X.append(IX(x,k));
#        P.append(IP(x,k));
                

# R = 300;
# XR = np.zeros((R,R));
# PR = np.zeros((R,R));
# for i in range(R):
#     for j in range(i,R):
#         XR[i,j] = X[j-i];
#         XR[j,i] = XR[i,j];
#         PR[i,j] = P[j-i];
#         PR[j,i] = PR[i,j];

# S=[];
# for n in range(1,R):
#     auxX = XR[0:n,0:n];
#     auxP = PR[0:n,0:n];
#     Cv = Msqrt(np.dot(auxX,auxP));
#     d1 = LA.eigvals(Cv+0.5*np.eye(n));
#     d2 = LA.eigvals(Cv-0.4999999999*np.eye(n));
#     S.append(np.dot(d1,np.log(d1))-np.dot(d2,np.log(d2)));

# plt.plot(S,'o-0')
# plt.show()

np.savez('500-sites', X, P);
