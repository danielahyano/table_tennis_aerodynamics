'''
This python program was made to compare graphs from the article https://pubs.aip.org/aapt/ajp/article/88/11/934/1058353/Flight-and-bounce-of-spinning-sports-balls Fig 11,
and mine.
'''
import numpy as np
import matplotlib.pyplot as plt


def readfile(filename):
    "Receives the file, and returns "
    plot_list = []
    with open(filename) as file:
        for item in file:
            x = float(item.split("\n")[0])
            plot_list.append(x)
    
    return plot_list

x_plot_noair = readfile("xNoAir.txt")
y_plot_noair = readfile("yNoAir.txt")
z_plot_noair = readfile("zNoAir.txt")

x_plot = readfile("xplot.txt")
y_plot = readfile("yplot.txt")
z_plot = readfile("zplot.txt")


plt.figure(0)
plt.clf()

plt.plot(x_plot, -1*np.array(y_plot),label="numerical, air resistance", color="blue")
plt.plot(x_plot_noair, -1*np.array(y_plot_noair), '--', label="analytic, vacuum", color='black')
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.xlim(0, 17)
plt.xticks(np.arange(0, 17.5, 2.5))
plt.ylim(-2.2, 2.2)
plt.yticks([-2,0,2], ["2", "0", "-2"])

plt.figure(1)
plt.clf()

plt.scatter(x_plot, z_plot, label="numerical, air resistance", color="blue")
plt.xlabel("x")
plt.ylabel("z")
plt.plot(x_plot_noair, z_plot_noair, '--',label="analytic, vacuum", color="black")
plt.xticks(np.arange(0, 20, 2.0))
plt.legend()
plt.xlim(0, 19)
plt.ylim(0, 2.2)


plt.show()