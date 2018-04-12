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

stats_data = clean_data(lines)
stats_data = stats_data[:len(t)]

stats_data2 = clean_data(lines2)
stats_data2 = stats_data2[:len(t)]


fig = plt.figure()
ax = fig.add_subplot(321)
ax.yaxis.grid(alpha=0.7, linestyle='dashed')
ax.xaxis.grid(alpha=0.7, linestyle='dashed')
plt.plot(t, autocorr(stats_data[:,0]),'b',label='Susceptible')
plt.legend()
plt.title('N = 100')
ax = fig.add_subplot(323)
ax.yaxis.grid(alpha=0.7, linestyle='dashed')
ax.xaxis.grid(alpha=0.7, linestyle='dashed')
plt.plot(t, autocorr(stats_data[:,1]),'y',label='Infected')
plt.legend()
ax = fig.add_subplot(325)
ax.yaxis.grid(alpha=0.7, linestyle='dashed')
ax.xaxis.grid(alpha=0.7, linestyle='dashed')
plt.plot(t, autocorr(stats_data[:,2]),'r', label='Recovered')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useOffset=False)
plt.legend()
ax = fig.add_subplot(322)
ax.yaxis.grid(alpha=0.7, linestyle='dashed')
ax.xaxis.grid(alpha=0.7, linestyle='dashed')
plt.plot(t, autocorr(stats_data2[:,0]),'b', label='Susceptible')
plt.legend()
plt.title('N = 500')
ax = fig.add_subplot(324)
ax.yaxis.grid(alpha=0.7, linestyle='dashed')
ax.xaxis.grid(alpha=0.7, linestyle='dashed')
plt.plot(t, autocorr(stats_data2[:,1]),'y', label='Infected')
plt.legend()
ax = fig.add_subplot(326)
ax.yaxis.grid(alpha=0.7, linestyle='dashed')
ax.xaxis.grid(alpha=0.7, linestyle='dashed')
plt.plot(t, autocorr(stats_data2[:,2]),'r', label='Recovered')
plt.legend()

plt.show()