import matplotlib.pyplot as plt
import numpy as np
import parameters as p
import atmosmodel as at
import init_functions as start

variables = ["w", "t", "qv", "qr", "qn"]
n = len(variables)

# Space conditions
z_initial = 0
z_final = 15
z_final = z_final / p.length_scale
space_steps = 101
dz = (z_final - z_initial) / (space_steps - 1)
space = [z_initial + i * dz for i in range(space_steps)]

# Time control
t_initial = 0
t_final = 2
t_final = t_final / p.time_scale

# Workspace
u = np.zeros((space_steps, n))
u[:, 0] = -1   # First Column Velocity
u[:, 1] = [start.initial_condition_theta(i, 0) for i in space]  # Second Column Temperature
u[:, 2] = [start.initial_condition_qv(i, 0) for i in space]  # Third Column QV
u[:, 3] = [start.initial_condition_qr(i, 0) for i in space]  # Fourth Column QR
u[:, 4] = [start.initial_condition_qn(i, 0) for i in space]  # Fiveth Column QN

# CFL criter
cfl = 0.9
dt = cfl * dz / np.max(np.abs(u[:, 0]))

basic_parameters = [p.g, u[0, 1], p.epsilon, p.B, u[0, 2]]
vt_parameters = [p.vt0, p.vtnd, p.q_star]

while t_initial < t_final:
    u = at.one_sided_model(dz, dt, u, space, basic_parameters, 0.75, vt_parameters)
    u[0, :] = u[1, :]  # Condiciones de Frontera Neumann
    u[space_steps - 1, :] = u[space_steps - 2, :]  # Condiciones de Frontera Neumann
    t_initial += dt

# bouyancy = [at.get_bouyancyforce(basic_parameters, i, u[0, 1], u[0, 2], u[0, 3]) for i in space]
# qvs = [at.approxfqv(i, u[0, 2]) for i in space]
# plt.plot(bouyancy)
# plt.plot(space, qvs)
print(np.shape(u))
plt.plot(u[:, 1], space)
plt.show()