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
z_plot_noair = readfile("zNoAir.txt")

x_plot = readfile("xplot.txt")
z_plot = readfile("zplot.txt")


v = np.array([readfile("vxplot.txt"), 
            readfile("vyplot.txt"), 
            readfile("vzplot.txt")])
v = np.swapaxes(v, 1, 0)
v1=v[0]
normv = round(np.linalg.norm(v1),1)

w = np.array([readfile("wxplot.txt"), 
            readfile("wyplot.txt"), 
            readfile("wzplot.txt")])
w = np.swapaxes(w, 1, 0)
w1 = w[0]
normw = np.linalg.norm(w1)
rps = round(normw/(2*np.pi),2)

backspin = False
if w1[1] < 0:
    backspin = True
plt.figure(0)
plt.clf()

# set table features
table_length = 2.7
table_height = 0.74
net_height = 0.152
# plot table
plt.plot(np.linspace(0, 2.7, 10), np.ones(10)*table_height, color="navy", linewidth=6)         #set table
plt.plot(np.ones(2)*(table_length/2), [table_height, table_height + net_height], color="green")  #set net
plt.plot(np.ones(2)*0.1, [0, table_height], color="blue")                                     #set first leg
plt.plot(np.ones(2)*(table_length-0.1), [0, table_height], color="blue")                      #set second leg

plt.scatter(x_plot, z_plot, c="crimson", s=5)
plt.xlabel("x (m)")
plt.ylabel("z (m)")

if backspin:
    plt.title(f'Backspin, initial velocity={normv} m/s, rps={rps}')
elif (normw == 0):
    plt.title(f"No spin, initial velocity={normv} m/s, rps={rps}")
else:
    plt.title(f"Topspin, initial velocity={normv} m/s, rps={rps}")

plt.xlim(-1, 3.7)
plt.ylim(0, 2)
plt.show()