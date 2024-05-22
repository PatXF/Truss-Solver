from truss import *
from nonlinear import nl_analysis

nodes = np.array([[0, 0], [4, 3], [8, 0]])
elements = np.array([[1, 2, 70, 645.2], [2, 3, 70, 645.2], [3, 1, 70, 645.2]])
boundary_conditions = {1: 0, 2: 0, 6: 0}
forces = {4: -2000}
cf, cd = solve(nodes, elements, boundary_conditions, forces)

print(cf)
print(cd)

final_disp, members_forces = nl_analysis(cd, cf, nodes, elements, 0.001)

print(final_disp)
print(members_forces)