import scipy as sc
import math as m
import matplotlib.pyplot as plt

beta = 100
gamma = 0.1
N = 500
s_init = 0.95
i_init = 0.05
t_init = 0.0

tmax = 150
dt = 0.01

t = sc.arange(0, tmax, dt)
s = sc.zeros(t.shape)
i = sc.zeros(t.shape)
r = sc.zeros(t.shape)

s[0] = s_init
i[0] = i_init

for k in range(1,len(t)):
    s[k] = s[k-1] - beta*i[k-1]*s[k-1]*dt / N
    i[k] = i[k-1] + (beta*i[k-1]*s[k-1]/N - gamma*i[k-1])*dt
    r[k] = r[k-1] + gamma*i[k-1]*dt


plt.plot(t,s,label='susceptible')
plt.plot(t,i,label='infected')
plt.plot(t,r,label='recovered')
plt.xlim([0,tmax])
plt.ylim([-0.1, 1])
plt.legend()
plt.show()

