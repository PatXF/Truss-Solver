import matplotlib.pyplot as plt
from helper import *


# draws an arrow and gives the magnitude : lab is the label for the arrow
def draw_arrow(x, y, mag, dirn, lab):
    mag /= 1000
    if mag == 0:
        return
    x1 = 0
    y1 = 0
    if dirn == "-x":
        x1 = -mag
    elif dirn == "x":
        x1 = mag
    elif dirn == "-y":
        y1 = -mag
    else:
        y1 = mag
    plt.annotate(lab, xy=(x + x1 + 0.1, y + y1 + 0.1))
    plt.arrow(x, y, x1, y1, color="blue", width=0.05, label=lab)


# provides a plot for the truss.
def visuals(nodes, forces, elements, ele_force):
    legends = []
    for i in range(len(elements)):
        x = []
        y = []
        x.append(nodes[int(elements[i][0] - 1)][0])
        x.append(nodes[int(elements[i][1] - 1)][0])
        y.append(nodes[int(elements[i][0] - 1)][1])
        y.append(nodes[int(elements[i][1] - 1)][1])
        plt.plot(x, y)
        try:
            legends.append(ele_force[i + 1])
        except KeyError:
            continue
    plt.title("Truss plot with joint forces as arrows and member forces as legends")
    plt.legend(legends)
    for key in forces.keys():
        if key % 2 == 1:
            node = int(key / 2)
            if forces[key] > 0:
                dirn = "x"
            else:
                dirn = "-x"
        else:
            node = int(key / 2) - 1
            if forces[key] > 0:
                dirn = "y"
            else:
                dirn = "-y"
        draw_arrow(nodes[node][0], nodes[node][1], abs(forces[key]), dirn, forces[key])
    plt.show()