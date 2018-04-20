import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('Final2b-CA.dat')

plt.imshow(data, cmap='binary')
plt.ylabel('Time')
plt.title('Rule 73 (Question 2b, 100 time steps)')
plt.show()