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

dt = 0.05
tmax = 2000

with open('reprod.dat', 'r') as d:
    lines = d.readlines()
d.close()

stats_data = clean_data(lines)
stats_data = stats_data[5:]
t = sc.arange(0, len(stats_data))

plt.plot(t,stats_data)
plt.show()