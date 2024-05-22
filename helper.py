import numpy as np
from visuals import *
from collections import OrderedDict


##Makes the code visually appealing by using arrows, legends etc.
def draw(nodes, complete_forces, elements, complete_disps, ele_forces):
    new_nodes = nodes.copy()
    for key in complete_forces.keys():
        if key % 2 == 0:
            index = int(key / 2) - 1
            new_nodes[index][1] += (complete_disps[key] / 50)
        else:
            index = int(key / 2)
            new_nodes[index][0] += (complete_disps[key] / 50)
    visuals(new_nodes, complete_forces, elements, ele_forces)


##Calculates the member forces, where the positive member force means compressive force and negative member force means
##tensile force. The calculations for member force is done as per the formulas for the non-linear geometric analysis.
def member_forces(element, nodes, disps):
    n1 = int(element[0])
    n2 = int(element[1])
    E = element[2]
    A = element[3]
    (x1, y1) = nodal_data(n1 - 1, nodes)
    (x2, y2) = nodal_data(n2 - 1, nodes)
    L = distance(x1, y1, x2, y2)
    x1_bar = x1 + disps[2 * n1 - 1]
    x2_bar = x2 + disps[2 * n2 - 1]
    y1_bar = y1 + disps[2 * n1]
    y2_bar = y2 + disps[2 * n2]
    L_bar = distance(x1_bar, y1_bar, x2_bar, y2_bar)
    Q = ((A * E) / L) * (L - L_bar)
    return Q


##Multiplies two matrices and returns a matrix
def mul(A, B):
    result = np.zeros((A.shape[0], B.shape[1]))
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            result[i][j] = (A[i, :] * B[:, j]).sum()
    return result


# Say the user didn't provide us with the zero forces at the free degrees of freedom and only provided us with the
# non-zero ones. This function finds all the possible forces at the free degrees of freedom.
def comp_forces(forces, boundary_conditions, num_nodes):
    for i in range(num_nodes):
        if i + 1 not in boundary_conditions.keys() and i + 1 not in forces.keys():
            forces[i + 1] = 0
    new_forces = OrderedDict(sorted(forces.items()))
    return new_forces


# When we represent S x D = F, when the known forces are not at the top of the F vector, it is known that transforming
# the rows of the Stiffness matrix and rows of the force vector would give us the known forces at the top. The code
# gives the top positions to the nearest known forces with a priority of FIFO.
def transform_row(stiffness, forces, num_nodes):
    ns = []
    order_f = []
    for key in forces.keys():
        order_f.append(key)
        ns.append(stiffness[key - 1])
    for i in range(num_nodes):
        if i + 1 not in forces.keys():
            order_f.append(i + 1)
            ns.append(stiffness[i])
    ns = np.array(ns)
    order_f = np.array(order_f)
    return ns, order_f


# When we represent S x D = F, when the unknown displacements are not at the top of the D vector, it is known that
# transforming the columns of the Stiffness matrix and rows of the displacement vector would give us the unknown
# displacements at the top. The code gives the top positions to the nearest unknown displacements with a priority of
# FIFO.
def transform_col(stiffness, boundary_conditions, num_nodes):
    ns = [[] for i in range(num_nodes)]
    order_d = []
    for i in range(num_nodes):
        if i + 1 not in boundary_conditions.keys():
            order_d.append(i + 1)
            t1 = stiffness[:, i]
            for j in range(num_nodes):
                ns[j].append(t1[j])
    for i in range(num_nodes):
        if i + 1 in boundary_conditions.keys():
            order_d.append(i + 1)
            t1 = stiffness[:, i]
            for j in range(num_nodes):
                ns[j].append(t1[j])
    ns = np.array(ns)
    order_d = np.array(order_d)
    return ns, order_d


# This function partitions the stiffness matrix into 4 parts
# 1st part corresponds to unknown displacements and known forces
# 2nd part corresponds to known displacements and known forces
# 3rd part corresponds to unknown displacements and unknown forces
# 4th part corresponds to known displacements and unknown forces
def partitions(stiff, fn, num):
    S1 = []
    for i in range(fn):
        temp = stiff[i, 0: fn]
        S1.append(temp)
    S2 = []
    for i in range(fn):
        temp = stiff[i, fn: num]
        S2.append(temp)
    S3 = []
    for i in range(fn, num):
        temp = stiff[i, 0: fn]
        S3.append(temp)
    S4 = []
    for i in range(fn, num):
        temp = stiff[i, fn: num]
        S4.append(temp)

    S1 = np.array(S1)
    S2 = np.array(S2)
    S3 = np.array(S3)
    S4 = np.array(S4)
    return S1, S2, S3, S4


# returns the x and y coordinate for a node
def nodal_data(node, nodes):
    return nodes[node][0], nodes[node][1]


# returns the euclidian distance between two points.
def distance(x1, y1, x2, y2):
    return np.sqrt(np.power((x1 - x2), 2) + np.power((y1 - y2), 2))