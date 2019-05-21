import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(-1,1,2000-1);
x = [];
y = [];
for i in range(len(t)):
	y.append(t[i]);
	x.append(2-np.sqrt(1-t[i]**2));

t = np.linspace(0,1,1000);
L = [];
for i in range(len(x)):
	x0 = x[i];
	y0 = y[i];
	aux1 = [];
	aux2 = [];
	for j in range(len(t)):
		aux1.append(np.sqrt( (y0*t[j])**2+ (1+t[j]*(-1+x0))**2 ) );
		aux2.append(np.arctan( t[j]*y0/( 1+t[j]*(-1+x0) ) ));
	res = 0;
	for j in range(len(aux1)-2):
		res = res+np.sqrt((aux1[j+1]-aux1[j])**2+aux1[j+1]**2*(aux2[j+1]-aux2[j])**2);
	L.append(res)

t = np.linspace(-1,1,2000-1);
plt.plot(t,L,'--')
plt.show()
