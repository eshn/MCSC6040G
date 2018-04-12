#### Eric Ng (100446517)
#### MSCS6040G Final Project
#### ODE.py - Plots of solution and simulation results

import scipy as sc
import math as m
import matplotlib.pyplot as plt

# Procedure to read and parse file contents
def clean_data(lines):
    data = []
    lines = [a.strip("\n,") for a in lines]
    for d in lines:
        data.append(d.split(" "))
    data = sc.double(sc.array(data))
    return data

N = 500

modelID = 6

# Parameters for cases
if modelID == 1: # Basic SIR
    beta = 0.55
    gamma = 0.1
    alpha = 0.0
    s_init = 0.9*N
    i_init = 0.1*N
    r_init = 0.0*N
    tmax = 50
    dt = 0.05
elif modelID == 2 or modelID == 3: # Basic SIRS (2 without waning population, 3 with waning population)
    beta = 0.55
    gamma = 0.1
    alpha = 0.1
    s_init = 0.9*N
    i_init = 0.1*N
    r_init = 0.0*N
    tmax = 50
    dt = 0.05
elif modelID == 4:
    tmax = 200
    dt = 0.05
elif modelID > 4: # SIRS with probability based recovery time and "bacteria like" movement
    tmax = 1000
    dt = 0.05

# Read from file
with open('SIRstat.dat', 'r') as d:
    lines = d.readlines()
d.close()

t = sc.arange(0, tmax, dt)
s = sc.zeros(t.shape)
i = sc.zeros(t.shape)
r = sc.zeros(t.shape)

stats_data = clean_data(lines)
stats_data = stats_data[:len(t)]

# Plots
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


    plt.plot(t,s,'b',label='Susceptible')
    plt.plot(t,i,'y',label='Infected')
    plt.plot(t,r,'r',label='Recovered')
    if modelID == 3:
        plt.plot(t,w,'k',label='Waning Population')


plt.plot(t,stats_data[:,0],'b--',label='Susceptible')
plt.plot(t,stats_data[:,1],'y--',label='Infected')
plt.plot(t,stats_data[:,2],'r--',label='Recovered')
# plt.plot(t,stats_data[:,0]+stats_data[:,2],'k--',label='S+R')
if modelID == 3:
    plt.plot(t,stats_data[:,3],'k--',label='W (Simulation)')
plt.xlim([0,tmax])
plt.legend()
plt.title('N = 100')
plt.xlabel('Time')
plt.ylabel('Population')
plt.show()

