import csv
import matplotlib.pyplot as plt

csvFiles = [
    'mpi_performance_2022_10_04_19_12_50.csv',
    'mpi_performance_2022_10_04_19_24_33.csv',
    'mpi_performance_2022_10_04_19_38_37.csv',
    'mpi_performance_2022_10_04_20_14_30.csv',
    'mpi_performance_2022_10_04_20_33_51.csv',
]

plotLines = {}

for filename in csvFiles:
    with open(filename, 'r') as csvfile:
        csvLines = csv.reader(csvfile, delimiter=',')
        next(csvLines)  # jump header line
        for row in csvLines:
            time_mpi = float(row[7])
            n = int(row[0])
            m = int(row[1])
            size = str(n) + 'x' + str(n)
            np = row[11]

            if np not in plotLines:
                plotLines[np] = {'x': [], 'y': []}
            try:
                i = plotLines[np]['x'].index(size)
                plotLines[np]['y'][i].append(time_mpi)
            except ValueError:
                plotLines[np]['x'].append(size)
                plotLines[np]['y'].append([time_mpi])

fig, ax = plt.subplots()
for key, value in plotLines.items():
    size = key
    x = value['x']

    y = value['y']
    y = [min(ys) for ys in y]   # min(ys), max(ys), numpy.average(ys)
    yerr = [
        [y[i]-min(ys) for i, ys in enumerate(value['y'])],
        [max(ys)-y[i] for i, ys in enumerate(value['y'])],
    ]
    ax.errorbar(x, y, yerr=yerr, marker='o', label=size, capsize=6, capthick=3)

plt.xlabel('Problem size')
plt.ylabel('Performance (in seconds)')
plt.title('MPI parallelization report', fontsize=20)
plt.grid()
plt.legend()
plt.show()

