import matplotlib.pyplot as plt  # Graficas
import numpy as np  # Vectores
import parameters as p  # Lista de parámetros
import atmosmodel as at  # Funciones principales
import init_functions as start  # Condiciones iniciales


def system_plots(velocity, temperature, vapor, water, core):
    wt_plots, ax = plt.subplots(1, 2)
    qvrn_plots, axq = plt.subplots(2, 2)
    ax[0].plot(velocity * p.velocity_scale, space * p.length_scale)
    ax[1].plot(temperature * p.temperature_scale, space * p.length_scale)
    axq[0, 0].plot(velocity * p.velocity_scale, space * p.length_scale)
    axq[0, 1].plot(vapor * p.ratio_scale, space * p.length_scale)
    axq[1, 0].plot(water * p.ratio_scale, space * p.length_scale)
    axq[1, 1].plot(core * p.ratio_scale, space * p.length_scale)
    plt.show()


def workspace_plots(wspace):
    system_plots(wspace[:, 0], wspace[:, 1], wspace[:, 2], wspace[:, 3], wspace[:, 4])


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
nz = 150  # Esta es nuesta m
space = np.linspace(z_initial, z_final, nz)
dz = space[1] - space[0]

# Time control
"""
Tiempo inicial y final
"""
t_initial = 0
t_initial = t_initial / p.time_scale
t_final = 10
t_final = t_final / p.time_scale

cut_omega = 3
cut_theta = 3
cut_qv = 3
cut_qr = 3
cut_qn = 3

cut_omega = cut_omega / p.length_scale  # No olvidemos adimensionalizar
cut_theta = cut_theta / p.temperature_scale
cut_qv = cut_qv / p.ratio_scale
cut_qr = cut_qr / p.ratio_scale
cut_qn = cut_qn / p.ratio_scale


# Workspace
v0 = 1
workspace = np.zeros((nz, 5))
workspace[:, 0] = v0
workspace[:, 1] = [start.heaviside(i - cut_theta) for i in space]
workspace[:, 2] = [start.heaviside(i - cut_qv) for i in space]
workspace[:, 3] = [start.heaviside(i - cut_qr) for i in space]
workspace[:, 4] = [start.heaviside(i - cut_qn) for i in space]

# CFL criter
cfl = 0.9
dt = cfl * dz / np.max(np.abs(workspace[:, 0]))

theta_0 = 300
theta_0 = theta_0 / p.temperature_scale

qv0 = 28
qv0 = qv0 / p.ratio_scale

basic_parameters = [p.g, theta_0, p.epsilon, p.B, qv0]
vt_parameters = [p.vt0, p.vtnd, p.q_star]

# workspace_plots(workspace)
result = at.resol_test(t_initial, t_final, cfl, dt, dz, workspace, vt_parameters, v0)
workspace_plots(result)

plt.show()
