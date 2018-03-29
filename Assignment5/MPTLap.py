# Computes MPT by Lap

import matplotlib.pyplot as plt
import numpy as np

def clean_data(lines, n):
    data = []
    lines = [a.strip("\n,") for a in lines]
    for i in range(len(lines)):
        if i % n != 0:
            data.append(lines[i].split(" "))
    data = np.array(data)
    return data

def readfile(filename):
    with open(filename, 'r') as d:
        data = d.readlines()
    d.close()
    return data

car = readfile('car.dat')
truck = readfile('truck.dat')

car = clean_data(car, 41)
truck = clean_data(truck, 6)

car = car[:,1]
truck = truck[:,1]

car = np.double(car)
truck = np.double(truck)

car = car.reshape((-1, 40))
truck = truck.reshape((-1, 5))

# columns are vehicle
# rows are position

car_time = np.zeros([20,40])
truck_time = np.zeros([10,5])

car_mean = np.zeros(20)
truck_mean = np.zeros(10)

for i in range(40):
    total = 0
    flag = 0 # flag to start recording to ensure for a full loop
    count = 0
    for j in range(car.shape[0]):
        if j == 0:
            continue
        elif flag == 0:
            if (abs(car[j,i] - car[j-1,i]) > 10):
                flag = 1
        else:
            if (abs(car[j,i] - car[j-1,i]) > 10):
                car_time[count,i] = total
                count += 1
                total = 0
            else:
                total += 1

for i in range(5):
    total = 0
    flag = 0
    count = 0
    for j in range(truck.shape[0]):
        if j == 0:
            continue
        elif flag == 0:
            if (abs(truck[j,i] - truck[j-1,i]) > 50):
                flag = 1
        else:
            if (abs(truck[j,i] - truck[j-1,i]) > 50):
                truck_time[count, i] += total
                count += 1
                total = 0
            else:
                total += 1

for i in range(len(truck_mean)):
    truck_mean[i] = np.mean(truck_time[i,:])

for i in range(len(car_mean)):
    car_mean[i] = np.mean(car_time[i,:])

truck_mean = truck_mean[:5]*0.02
car_mean = car_mean[:15]*0.02

x = np.linspace(1,15,15)
xs = np.linspace(0,15)
xs2 = np.linspace(0,5)
plt.subplot(121)
plt.scatter(x, car_mean, label='Mean Passage Time by Lap ')
plt.plot(xs, 5.766*np.ones(xs.shape),'r--', label='Overall Mean Passage Time')
plt.xlim([0, 15])
plt.ylim([4, 8])
plt.xlabel('Mean Passage Time')
plt.ylabel('Lap')
plt.title('Mean Passage Time by Lap (Small Vehicles)')
plt.legend()
plt.subplot(122)
plt.scatter(x[:5], truck_mean, label='Mean Passage Time by Lap ')
plt.plot(xs2, 18.483*np.ones(xs2.shape),'r--', label='Overall Mean Passage Time')
plt.xlim([0, 5])
plt.ylim([17, 20])
plt.xlabel('Mean Passage Time')
plt.ylabel('Lap')
plt.title('Mean Passage Time by Lap (Large Vehicles)')
plt.legend()
plt.show()

print(truck_mean)
print(car_mean)