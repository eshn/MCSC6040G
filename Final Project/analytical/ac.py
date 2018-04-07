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

def autocorr(x):
    result = sc.correlate(x, x, mode='full')
    result = result[len(x)-1:]
    return result

with open('model6.dat', 'r') as d:
    lines = d.readlines()
d.close()

with open('model7.dat', 'r') as d:
    lines2 = d.readlines()
d.close()

t = sc.arange(0, 1000, 0.05)
s = sc.zeros(t.shape)
i = sc.zeros(t.shape)
r = sc.zeros(t.shape)

stats_data = clean_data(lines)
stats_data = stats_data[:len(t)]

stats_data2 = clean_data(lines2)
stats_data2 = stats_data2[:len(t)]

fig = plt.figure()
ax = fig.add_subplot(231)
ax.yaxis.grid(alpha=0.7, linestyle='dashed')
ax.xaxis.grid(alpha=0.7, linestyle='dashed')
plt.plot(t, autocorr(stats_data[:,0]),'b',label='Susceptible')
plt.legend()
plt.ylabel('N = 100')
ax = fig.add_subplot(232)
ax.yaxis.grid(alpha=0.7, linestyle='dashed')
ax.xaxis.grid(alpha=0.7, linestyle='dashed')
plt.plot(t, autocorr(stats_data[:,1]),'y',label='Infected')
plt.legend()
ax = fig.add_subplot(233)
ax.yaxis.grid(alpha=0.7, linestyle='dashed')
ax.xaxis.grid(alpha=0.7, linestyle='dashed')
plt.plot(t, autocorr(stats_data[:,2]),'r', label='Recovered')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useOffset=False)
plt.legend()
ax = fig.add_subplot(234)
ax.yaxis.grid(alpha=0.7, linestyle='dashed')
ax.xaxis.grid(alpha=0.7, linestyle='dashed')
plt.plot(t, autocorr(stats_data2[:,0]),'b', label='Susceptible')
plt.legend()
plt.ylabel('N = 500')
ax = fig.add_subplot(235)
ax.yaxis.grid(alpha=0.7, linestyle='dashed')
ax.xaxis.grid(alpha=0.7, linestyle='dashed')
plt.plot(t, autocorr(stats_data2[:,1]),'y', label='Infected')
plt.legend()
ax = fig.add_subplot(236)
ax.yaxis.grid(alpha=0.7, linestyle='dashed')
ax.xaxis.grid(alpha=0.7, linestyle='dashed')
plt.plot(t, autocorr(stats_data2[:,2]),'r', label='Recovered')
plt.legend()

plt.show()