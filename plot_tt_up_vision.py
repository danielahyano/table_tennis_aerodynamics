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
v = np.array([readfile("vxplot.txt"), 
            readfile("vyplot.txt"), 
            readfile("vzplot.txt")])
v = np.swapaxes(v, 1, 0)
fig, ax = plt.subplots()
scat = ax.scatter(x_plot[0], y_plot[0], c="red", s=5, label=f"|v|: {round(np.linalg.norm(v[0]), 2)} m/s")

ax.set(xlim=[-1, 3.7], ylim=[-0.5, 2])
L=plt.legend(loc=1) #Define legend objects

def update(frame):
    # for each frame, update the data stored on each artist.
    x = x_plot[:frame]
    y = y_plot[:frame]
    lab = f"|v|: {round(np.linalg.norm(v[len(x)]), 2)} m/s"
    # update the scatter plot:
    data = np.stack([x, y]).T
    scat.set_offsets(data)
    L.get_texts()[0].set_text(lab) #Update label each at frame
    return scat

# set table features
table_length = 2.7
table_width = 1.525
table_height = 0.74
net_height = 0.152
# plot table
plt.plot(np.linspace(0, table_length, 10), np.ones(10)*table_width, color="blue")       
plt.plot(np.linspace(0, table_length, 10), np.zeros(10), color="blue")   
plt.plot(np.linspace(0, table_length, 10), np.ones(10)*table_width/2, color="darkviolet")         
plt.plot(np.zeros(10), np.linspace(0, table_width, 10), color="blue")  
plt.plot(np.ones(10)*table_length, np.linspace(0, table_width, 10), color="blue")  
plt.plot(np.ones(10)*table_length/2, np.linspace(0, table_width, 10), color="green")  
ani = animation.FuncAnimation(fig=fig, func=update, frames=len(x_plot), interval=1)
plt.show()