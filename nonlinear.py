from helper import *
from gj import *


# calculated the change from the previous iteration and the delta displacement obtained
def calc_change(d, dd):
    num = 0
    den = 0
    for i in range(len(d)):
        num += (dd[i] ** 2)
        den += (d[i][0] ** 2)
    return np.sqrt(num/den)


# similar to linear analysis but contains the displacement of the nodes are large enough to take part in calculations.
def nl_analysis(disps, forces, nodes, elements, tolerance):
    nodes_to_look_for = []
    d = []
    P = []
    for key in disps.keys():
        if disps[key] != 0:
            nodes_to_look_for.append(key)
            d.append([disps[key]])
            P.append([forces[key]])
    f = np.zeros((len(nodes_to_look_for), 1))
    num_nodes = nodes.shape[0]
    global_tangent_stiffness = np.zeros((2 * num_nodes, 2 * num_nodes))
    for element in elements:
        n1 = int(element[0])
        n2 = int(element[1])
        E = element[2]
        A = element[3]
        (x1, y1) = nodal_data(n1 - 1, nodes)
        (x2, y2) = nodal_data(n2 - 1, nodes)
        x1_bar = x1 + disps[2*n1 - 1]
        x2_bar = x2 + disps[2*n2 - 1]
        y1_bar = y1 + disps[2*n1]
        y2_bar = y2 + disps[2*n2]
        L_bar = distance(x1_bar, y1_bar, x2_bar, y2_bar)
        cx = (x2_bar - x1_bar) / L_bar
        cy = (y2_bar - y1_bar) / L_bar
        Kt = np.array([[cx ** 2, cx * cy, -(cx ** 2), -cx * cy],
                       [cx * cy, cy ** 2, -cx * cy, -(cy ** 2)],
                       [-(cx ** 2), -cx * cy, cx ** 2, cx * cy],
                       [-cx * cy, -(cy ** 2), cx * cy, cy ** 2]])
        L = distance(x1, y1, x2, y2)
        Kt = ((E * A) / L) * Kt
        g = np.array([[-(cy ** 2), cx * cy, cy ** 2, -cx * cy],
                      [cx * cy, -(cx ** 2), -cx * cy, cx ** 2],
                      [cy ** 2, -cx * cy, -(cy ** 2), cx * cy],
                      [-cx * cy, cx ** 2, cx * cy, -(cx ** 2)]])
        Q = ((A * E) / L) * (L - L_bar)
        Kt = Kt + (Q * g) / L_bar
        l1 = 2 * n1 - 1
        l2 = 2 * n1
        l3 = 2 * n2 - 1
        l4 = 2 * n2
        indexes = [l1, l2, l3, l4]
        for i in range(4):
            for j in range(4):
                global_tangent_stiffness[indexes[i] - 1][indexes[j] - 1] += Kt[i][j]
        F = np.array([cx, cy, -cx, -cy]) * Q
        for i in range(len(nodes_to_look_for)):
            for j in range(4):
                if nodes_to_look_for[i] == indexes[j]:
                    f[i][0] += F[j]
                    break

    S = []
    for i in range(len(nodes_to_look_for)):
        temp = []
        for j in range(len(nodes_to_look_for)):
            temp.append(global_tangent_stiffness[nodes_to_look_for[i] - 1][nodes_to_look_for[j] - 1])
        S.append(temp)
    S = np.array(S)
    dU = P - f
    dd = gj(S, dU)
    if calc_change(d, dd) <= tolerance:
        c = 0
        for key in disps.keys():
            if disps[key] != 0:
                disps[key] += dd[c]
                c += 1
        mem_for = {}
        c = 1
        for element in elements:
            q = member_forces(element, nodes, disps)
            mem_for[c] = q
            c += 1
        draw(nodes, forces, elements, disps, mem_for)
        return disps, mem_for

    c = 0
    for key in disps.keys():
        if disps[key] != 0:
            disps[key] += dd[c]
            c += 1
    return nl_analysis(disps, forces, nodes, elements, tolerance)