'''
Author: Daniela Yano
Date: April 2024

This code is designed to plot the results from 

'''


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

wx_plot = readfile("wxplot.txt")
wy_plot = readfile("wyplot.txt")
wz_plot = readfile("wzplot.txt")

v = np.array([readfile("vxplot.txt"), 
            readfile("vyplot.txt"), 
            readfile("vzplot.txt")])
v = np.swapaxes(v, 1, 0)

fig, ax = plt.subplots()

scat = ax.scatter(x_plot[0], z_plot[0], c="pink", s=5, label=f"|v|: {round(np.linalg.norm(v[0]), 2)} m/s, wx = {wx_plot[0]}, wy = {wy_plot[0]}, wz = {wz_plot[0]}")
ax.set(xlim=[-1, 3.7], ylim=[0, 2])
L=plt.legend(loc=1) #Define legend objects
# ax.legend()

def update(frame):
    # for each frame, update the data stored on each artist.
    x = x_plot[:frame]
    z = z_plot[:frame]
    lab = f"|v|: {round(np.linalg.norm(v[len(x)]), 2)} m/s, , wx = {wx_plot[len(x)]}, wy = {wy_plot[len(x)]}, wz = {wz_plot[len(x)]}"

    # update the scatter plot:
    data = np.stack([x, z]).T
    scat.set_offsets(data)
    L.get_texts()[0].set_text(lab) #Update label each at frame
    return scat

# set table features
table_length = 2.7
table_height = 0.74
net_height = 0.152
# plot table
plt.plot(np.linspace(0, 2.7, 10), np.ones(10)*table_height, color="blue")                             #set table
plt.plot(np.ones(2)*(table_length/2), [table_height, table_height + net_height], color="green")  #set net
plt.plot(np.ones(2)*0.1, [0, table_height], color="blue")                                     #set first leg
plt.plot(np.ones(2)*(table_length-0.1), [0, table_height], color="blue")                      #set second leg

ani = animation.FuncAnimation(fig=fig, func=update, frames=len(x_plot), interval=1)
plt.show()