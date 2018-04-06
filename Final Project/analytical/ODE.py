import scipy as sc
import math as m
import matplotlib.pyplot as plt

def clean_data(lines):
    data = []
    lines = [a.strip("\n,") for a in lines]
    for d in lines:
        data.append(d.split(" "))
    data = sc.double(sc.array(data))
    return data

N = 500
modelID = 5

if modelID == 1: # Basic SIR
    beta = 0.6
    gamma = 0.12
    alpha = 0.0
    s_init = 0.95*N
    i_init = 0.05*N
    r_init = 0.0*N
    tmax = 50
    dt = 0.05
elif modelID == 2 or modelID == 3: # Basic SIRS (2 without waning population, 3 with waning population)
    beta = 0.6
    gamma = 0.13
    alpha = 0.12
    s_init = 0.9*N
    i_init = 0.1*N
    r_init = 0.0*N
    tmax = 50
    dt = 0.05
elif modelID == 4: # SIRS with probability based recovery time and "bacteria like" movement
    tmax = 2000
    dt = 0.05
elif modelID == 5:  # SIRS with probability based recovery time and "bacteria like" movement
    tmax = 2000
    dt = 0.05

with open('SIRstat.dat', 'r') as d:
    lines = d.readlines()
d.close()

t = sc.arange(0, tmax, dt)
s = sc.zeros(t.shape)
i = sc.zeros(t.shape)
r = sc.zeros(t.shape)

stats_data = clean_data(lines)
stats_data = stats_data[:len(t)]

if modelID < 4:
    s[0] = s_init
    i[0] = i_init
    if modelID == 3:
        w = sc.zeros(t.shape)
        for k in range(1,len(t)):
            s[k] = s[k-1] - beta*i[k-1]*s[k-1]*dt / N
            i[k] = i[k-1] + (beta*i[k-1]*(s[k-1]+w[k-1])/N - gamma*i[k-1])*dt
            r[k] = r[k-1] + gamma*i[k-1]*dt - alpha*r[k-1]*dt
            w[k] = w[k-1] - beta*w[k-1]*i[k-1]*dt / N + alpha*r[k-1]*dt
    else:
        for k in range(1,len(t)):
            s[k] = s[k-1] - beta*i[k-1]*s[k-1]*dt / N + alpha*r[k-1]*dt
            i[k] = i[k-1] + (beta*i[k-1]*s[k-1]/N - gamma*i[k-1])*dt
            r[k] = r[k-1] + gamma*i[k-1]*dt - alpha*r[k-1]*dt


    plt.plot(t,s,'b',label='S')
    plt.plot(t,i,'y',label='I')
    plt.plot(t,r,'r',label='R')
    if modelID == 3:
        plt.plot(t,w,'k',label='W')


plt.plot(t,stats_data[:,0],'b-.',label='Ssim')
plt.plot(t,stats_data[:,1],'y-.',label='Isim')
plt.plot(t,stats_data[:,2],'r-.',label='Rsim')
if modelID == 3:
    plt.plot(t,stats_data[:,3],'k-.',label='Wsim')
plt.xlim([0,tmax])
# plt.ylim([-0.1, 1])
plt.legend()
plt.show()

