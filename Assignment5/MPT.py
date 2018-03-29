# Computes Mean Passage Time

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

car_time = np.zeros(car.shape[1])
truck_time = np.zeros(truck.shape[1])

car_loop_time = np.zeros(20)

for i in range(40):
    total = 0
    flag = 0
    count = 0
    for j in range(car.shape[0]):
        if j == 0:
            continue
        elif flag == 0:
            if (abs(car[j,i] - car[j-1,i]) > 10):
                flag = 1
        else:
            if (abs(car[j,i] - car[j-1,i]) > 10):
                count += 1
                car_time[i] += total
                total = 0
            else:
                total += 1
    car_time[i] = car_time[i] / count

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
                count += 1
                truck_time[i] += total
                total = 0
            else:
                total += 1

    truck_time[i] = truck_time[i] / count

print(np.mean(car_time)*0.02)
print(np.mean(truck_time)*0.02)