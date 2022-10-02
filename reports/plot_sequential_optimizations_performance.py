import matplotlib.pyplot as plt
import csv

plotLines = {}

with open('sequential_optimizations_performance_2022_10_02_13_59_31.csv', 'r') as csvfile:
    csvLines = csv.reader(csvfile, delimiter=',')
    next(csvLines)
    for row in csvLines:
        time_mpi = float(row[7])
        np = int(row[11])
        opt = row[12]

        if opt not in plotLines:
            plotLines[opt] = {'x': [], 'y': []}
        plotLines[opt]['x'].append(np)
        plotLines[opt]['y'].append(time_mpi)

for opt, xy in plotLines.items():
    opt = "-0" + opt
    plt.plot(xy['x'], xy['y'], marker='o', label=opt)

plt.xlabel('Number of busy cores in a single node')
plt.ylabel('Performance (in seconds)')
plt.title('Optimizations performance report', fontsize=20)
plt.grid()
plt.legend()
plt.show()

