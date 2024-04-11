import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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
z_plot = readfile("zplot.txt")


fig, ax = plt.subplots()

scat = ax.scatter(x_plot[0], y_plot[0], c="pink", s=5)
ax.set(xlim=[-1, 3.7], ylim=[-1, 2], xlabel='X [m]', ylabel='Y [m]')
# ax.legend()

def update(frame):
    # for each frame, update the data stored on each artist.
    x = x_plot[:frame]
    y = y_plot[:frame]
    # update the scatter plot:
    data = np.stack([x, y]).T
    scat.set_offsets(data)
    return scat

# set table features
table_length = 2.7
table_width = 1.525
table_height = 0.74
net_height = 0.152
# plot table
plt.plot(np.linspace(0, table_length, 10), np.ones(10)*table_width, color="blue")       
plt.plot(np.linspace(0, table_length, 10), np.zeros(10), color="blue")   
plt.plot(np.linspace(0, table_length, 10), np.ones(10)*table_width/2, color="blue")         
plt.plot(np.zeros(10), np.linspace(0, table_width, 10), color="blue")  
plt.plot(np.ones(10)*table_length, np.linspace(0, table_width, 10), color="blue")  
plt.plot(np.ones(10)*table_length/2, np.linspace(0, table_width, 10), color="blue")  
ani = animation.FuncAnimation(fig=fig, func=update, frames=len(x_plot), interval=1)
plt.show()