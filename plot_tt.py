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

x_plot = readfile("xplot.txt")
y_plot = readfile("yplot.txt")


table_length = 2.7
table_height = 0.76
net_height = 0.152

plt.figure(0)
plt.clf()

# plot table
plt.plot(np.linspace(0, 2.7, 10), np.ones(10)*0.76, color="blue")                             #set table
plt.plot(np.ones(2)*(table_length/2), [table_height, table_height + net_height], color="green")  #set net
plt.plot(np.ones(2)*0.1, [0, table_height], color="blue")                                     #set first leg
plt.plot(np.ones(2)*(table_length-0.1), [0, table_height], color="blue")                      #set second leg

# plt.plot(x_plot_noair, y_plot_noair)
# plt.plot(x_plot, y_plot)
plt.xlabel("x")
plt.ylabel("y")

plt.xlim(-1, 3.7)
plt.ylim(0, 3)
plt.show()