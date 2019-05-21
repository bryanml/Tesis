import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(-1,1,2000);
L = np.zeros(len(t));
for i in range(len(t)):
	L[i] = np.sqrt(0.25*(5-np.sqrt(1-(t[i]**2))));

plt.plot(t,L,'--')
plt.show()
