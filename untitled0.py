import numpy as np
import matplotlib.pyplot as plt


print('Amplitud de ángulo - Fracción de pi:')

f = input();

print("cut-off rc:")

rc =input();
rc = float(rc);

h = 0.0001;
arc = float(f)*np.pi;
do = arc/2; ### Dtheta /2, media longitud de arco
c =np.sqrt(((np.sinh(rc)**2)/np.tan(do)**2)/(np.cosh(rc)**2+1/np.tan(do)**2)); #### sinh(r0)
r0 = np.arcsinh(c);
N = int((rc-r0)/h);
o = np.zeros(N);
r = [r0];
o[0] = 0;
for i in range(1,N):
	o[i] = o[i-1]+h*c/(np.sinh(r[i-1]+h)*np.sqrt(np.sinh(r[i-1]+h)**2-c**2));
	r.append(r0+i*h);
    
o2 = np.zeros(len(o));
o3 = np.zeros(len(o));
o4 = np.zeros(len(o));
for j in range(len(o)):
	o2[j] = -o[j];
	o3[j] = o[j]+np.pi;
	o4[j] = o2[j]+np.pi;

plt.polar(o,r,'--k',o2,r,'--k',o3,r,'--k',o4,r,'--k')
plt.show()

L1 = [];
L2 = [];
t = np.linspace(0,1,1000);
for i in range(len(o)):
    c1 = r[i]*np.cos(o[i]);
    c2 = r[i]*np.sin(o[i]);
    aux1 = np.zeros(len(t));
    aux2 = np.zeros(len(t));
    a1 = r[i]*np.cos(o2[i]);
    a2 = r[i]*np.sin(o2[i]);
    aux3 = np.zeros(len(t));
    aux4 = np.zeros(len(t));
    for j in range(len(t)):
        aux1[j] = np.sqrt((t[j]*(c1+r0)-r0)**2+t[j]**2*c2**2);
        aux2[j] = np.arctan(t[j]*c2/(t[j]*(c1+r0)-r0));
        aux3[j] = np.sqrt((t[j]*(a1+r0)-r0)**2+t[j]**2*a2**2);
        aux4[j] = np.arctan(t[j]*a2/(t[j]*(a1+r0)-r0));  
    res1 = 0;
    res2 = 0;
    for j in range(len(aux1)-2):
        res1 = res1+np.sqrt((aux1[j+1]-aux1[j])**2+np.sinh(aux1[j+1])**2*(aux2[j+1]-aux2[j])**2);
        res2 = res2+np.sqrt((aux3[j+1]-aux3[j])**2+np.sinh(aux3[j+1])**2*(aux4[j+1]-aux4[j])**2);
    L2.append(res2);
    L1.append(res1);

plt.plot(o,L1,'--k',o2,L2,'--k')
plt.show()
