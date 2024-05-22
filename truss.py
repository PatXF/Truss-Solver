from helper import *
from gj import *
from visuals import *
import numpy as np

member_stiffness = np.array([[1, 0, -1, 0], [0, 0, 0, 0],
                             [-1, 0, 1, 0], [0, 0, 0, 0]])


# solves the truss using partitions and linear analysis.
def truss_solver(nodes, elements, forces, boundary_conditions):
    num_nodes = nodes.shape[0]
    global_stiffness = np.zeros((2*num_nodes, 2*num_nodes))
    for element in elements:
        n1 = int(element[0])
        n2 = int(element[1])
        E = element[2]
        A = element[3]
        (x1, y1) = nodal_data(n1 - 1, nodes)
        (x2, y2) = nodal_data(n2 - 1, nodes)
        L = distance(x1, y1, x2, y2)
        cos_theta = (x2 - x1) / L
        sin_theta = (y2 - y1) / L
        ms = ((E * A) / L) * member_stiffness
        T = np.array([[cos_theta, -sin_theta, 0, 0],
                      [sin_theta, cos_theta, 0, 0],
                      [0, 0, cos_theta, -sin_theta],
                      [0, 0, sin_theta, cos_theta]])
        ms = mul(mul(T, ms), T.T)
        l1 = 2 * n1 - 1
        l2 = 2 * n1
        l3 = 2 * n2 - 1
        l4 = 2 * n2
        indexes = [l1, l2, l3, l4]
        for i in range(4):
            for j in range(4):
                global_stiffness[indexes[i] - 1][indexes[j] - 1] += ms[i][j]
    global_stiffness, new_f_order = transform_row(global_stiffness, forces, 2*num_nodes)
    global_stiffness, new_d_order = transform_col(global_stiffness, boundary_conditions, 2*num_nodes)
    S1, S2, S3, S4 = partitions(global_stiffness, len(forces), 2*num_nodes)
    known = [[boundary_conditions[key]] for key in boundary_conditions.keys()]
    known = np.array(known)
    prod = mul(S2, known)
    known_forces = []
    for key in forces.keys():
        known_forces.append([forces[key]])
    known_forces = np.array(known_forces)
    for i in range(len(known_forces)):
        known_forces[i][0] -= prod[i][0]
    unknown_disp = gj(S1, known_forces)
    ud = []
    for i in range(len(unknown_disp)):
        ud.append([unknown_disp[i]])
    ud = np.array(ud)
    ans1 = mul(S3, ud)
    ans2 = mul(S4, known)
    unknown_forces = []
    for i in range(len(ans1)):
        unknown_forces.append([ans1[i][0] + ans2[i][0]])

    return unknown_disp, unknown_forces, new_d_order, new_f_order


# method to solve the truss and present a visually appealing result
def solve(nodes, elements, boundary_conditions, forces):
    temp_vis = {}
    visuals(nodes, forces, elements, temp_vis)
    forces = comp_forces(forces, boundary_conditions, nodes.shape[0] * 2)
    complete_forces = {}
    complete_disps = {}
    u_disp, u_force, order_disp, order_force = truss_solver(nodes, elements, forces, boundary_conditions)
    for key in boundary_conditions.keys():
        complete_disps[key] = boundary_conditions[key]
    for key in forces.keys():
        complete_forces[key] = forces[key]
    known_disp_index_len = len(boundary_conditions)
    known_forces_index_len = len(forces)
    for i in range(2 * len(nodes) - known_disp_index_len):
        complete_disps[order_disp[i]] = u_disp[i]
    for i in range(known_forces_index_len, 2 * len(nodes)):
        complete_forces[order_force[i]] = u_force[i - known_forces_index_len][0]

    mem_for = {}
    c = 1
    for element in elements:
        q = member_forces(element, nodes, complete_disps)
        mem_for[c] = q
        c += 1

    draw(nodes, complete_forces, elements, complete_disps, mem_for)
    return complete_forces, complete_disps