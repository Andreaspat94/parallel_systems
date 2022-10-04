import csv
import matplotlib.pyplot as plt
import numpy

csvFiles = [
    'sequential_optimizations_performance_2022_10_02_13_59_31.csv',
    'sequential_optimizations_performance_2022_10_02_15_47_09.csv',
]

plotLines = {}

for filename in csvFiles:
    with open(filename, 'r') as csvfile:
        csvLines = csv.reader(csvfile, delimiter=',')
        next(csvLines)
        for row in csvLines:
            time_mpi = float(row[7])
            np = int(row[11])
            opt = row[12]

            if opt not in plotLines:
                plotLines[opt] = {'x': [], 'y': []}
            try:
                i = plotLines[opt]['x'].index(np)
                plotLines[opt]['y'][i].append(time_mpi)
            except ValueError:
                plotLines[opt]['x'].append(np)
                plotLines[opt]['y'].append([time_mpi])

fig, ax = plt.subplots()
for key, value in plotLines.items():
    opt = "-0" + key
    x = value['x']

    y = value['y']
    y = [min(ys) for ys in y]   # min(ys), max(ys), numpy.average(ys)
    yerr = [
        [y[i]-min(ys) for i, ys in enumerate(value['y'])],
        [max(ys)-y[i] for i, ys in enumerate(value['y'])],
    ]
    ax.errorbar(x, y, yerr=yerr, marker='o', label=opt, capsize=6, capthick=3)

plt.xlabel('Number of busy cores in a single node')
plt.ylabel('Performance (in seconds)')
plt.title('Optimizations performance report', fontsize=20)
plt.grid()
plt.legend()
plt.show()

