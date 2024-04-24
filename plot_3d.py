import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

def readfile(filename):
    "Receives the file, and returns "
    plot_list = []
    with open(filename) as file:
        for item in file:
            x = float(item.split("\n")[0])
            plot_list.append(x)
    
    return np.array(plot_list)
m=-300
x_data = readfile("xplot.txt")[:m]
y_data = -1*readfile("yplot.txt")[:m]
z_data = readfile("zplot.txt")[:m]
num_points = np.shape(x_data)[0]

# Create a figure and a 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Setting the axes properties
ax.set(xlim3d=(-1, 3))
ax.set(ylim3d=(-1, 3))
ax.set(zlim3d=(0, 1.3))
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

# set table features
table_length = 2.7
table_height = 0.74
net_height = 0.152

# Define width and length of the table
width = 1.525
length = 2.7
height = 0.74

# Define data for the table
x = np.arange(0, width, 0.1)
y = np.arange(0, length, 0.1)
z = np.ones((len(y), len(x)))*height

# Create meshgrid for 3D plotting
X, Y = np.meshgrid(x, y)

# Plot table surface
ax.plot_surface(Y, X, z, alpha=0.5)
leg_height = np.linspace(0, height, 2)  # Two points for the height of each leg
offset = 0.1
leg_positions = [(0+offset, 0+offset), (0+offset, width-offset), (length-offset-0.1, 0+offset), (length-offset-0.1, width-offset)]  # Positions of the legs
for pos in leg_positions:
    ax.plot([pos[0], pos[0]], [pos[1], pos[1]], leg_height, color='k')  # Plotting legs

# Plot net in the middle
net_x = width / 2
net_y = length / 2
net_width = width
net_height = 0.125
ax.plot([net_y, net_y], [net_x - net_width / 2, net_x + net_width / 2], [height+net_height, height+net_height], color='b')
ax.plot([net_y, net_y], [net_y, net_y], np.linspace(height, height+net_height, 2), color='black' )
ax.plot([net_y, net_y], [0.1, 0.1], np.linspace(height, height+net_height, 2), color='black' )
ax.scatter(x_data, y_data, z_data, c='crimson', marker='o')
ax.legend()

plt.show()

import pickle
pickle.dump(fig, open('FigureObject.fig.pickle', 'wb'))