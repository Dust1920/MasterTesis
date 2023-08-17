import matplotlib.pyplot as plt  # Graficas
import numpy as np  # Vectores
import parameters as p  # Lista de parámetros
import atmosmodel as at  # Funciones principales
import init_functions as start  # Condiciones iniciales

variables = ["w", "t", "qv", "qr", "qn"]
n = len(variables)

# Space conditions
"""
Altura inicial y final
"""
z_initial = 0
z_initial = z_initial / p.length_scale
z_final = 15
z_final = z_final / p.length_scale

"""
Numero de pasos, delta z y vector de espacio
"""
space_steps = 101
dz = (z_final - z_initial) / (space_steps - 1)
space = [z_initial + i * dz for i in range(space_steps)]

# Time control
"""
Tiempo inicial y final
"""
t_initial = 0
t_initial = t_initial / p.time_scale
t_final = 0.2
t_final = t_final / p.time_scale


# Workspace
u = np.zeros((space_steps, n))
u[:, 0] = -1   # First Column Velocity
u[:, 1] = [start.initial_function(0, 1, start.s_theta, i, 0) for i in space]  # Second Column Temperature
u[:, 2] = [start.initial_function(0, 1, start.s_qv, i, 0) for i in space]  # Third Column QV
u[:, 3] = [start.initial_function(0, 1, start.s_qr, i, 0) for i in space]  # Fourth Column QR
u[:, 4] = [start.initial_function(0, 1, start.s_qn, i, 0) for i in space]  # Fiveth Column QN

# CFL criter
cfl = 0.9
dt = cfl * dz / np.max(np.abs(u[:, 0]))

basic_parameters = [p.g, u[0, 1], p.epsilon, p.B, u[0, 2]]
vt_parameters = [p.vt0, p.vtnd, p.q_star]

while t_initial < t_final:
    u = at.one_sided_model_2(dz, dt, u, vt_parameters, p.B, u[0, 1], u[0, 2])
    u[0, :] = u[1, :]  # Condiciones de Frontera Neumann
    u[space_steps - 1, :] = u[space_steps - 2, :]  # Condiciones de Frontera Neumann
    dt = cfl * dz / np.max(np.abs(u[:, 0]))
    t_initial += dt

"""
while t_initial < t_final:
    u = at.one_sided_model(dz, dt, u, space, basic_parameters, 0.75, vt_parameters)
    u[0, :] = u[1, :]  # Condiciones de Frontera Neumann
    u[space_steps - 1, :] = u[space_steps - 2, :]  # Condiciones de Frontera Neumann
    t_initial += dt
"""

# bouyancy = [at.get_bouyancyforce(basic_parameters, i, u[0, 1], u[0, 2], u[0, 3]) for i in space]
# qvs = [at.approxfqv(i, u[0, 2]) for i in space]
# plt.plot(bouyancy)
# plt.plot(space, qvs)
print(np.shape(u))
plt.plot(u[:, 1], space)
plt.show()
