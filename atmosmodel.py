import numpy as np


def approxfqv(z, c):
    # Dr.Gerardo Hernandez Dueñas
    a = 18.04
    b = 3.27
    c0 = 0.1
    d = 3.48
    pz = (1 - b * np.log(1 + c0 * z)) ** d
    fqv = (c / pz) * np.exp(- a * (1 / ((1 - b * np.log(1 + c0 * z)) * (1 + c0 * z)) - 1))
    return fqv


def get_bouyancyforce(par, z, theta, qv, qr):
    g = par[0]
    theta_0 = par[1]
    epsilon = par[2]
    b0 = par[3]
    qv_0 = par[4]
    theta_hat = theta_0 + b0 * z
    qv_hat = approxfqv(z, qv_0)
    b = g / theta_0 * ((theta - theta_hat) + epsilon * theta_0 * (qv - qv_hat) / 100 - theta_0 * qr / 100)
    return b


def get_scale_time_condesation(tau_0, qn0, qn, gamma):
    tau_c = tau_0 * np.exp(((qn - qn0) / gamma) ** 2)
    return tau_c


def get_condensation(tau_0, qn0, qn, gamma, qv, qvs0, z):
    cd = 1 / get_scale_time_condesation(tau_0, qn0, qn, gamma) * np.max([qv - approxfqv(z, qvs0), 0])
    return cd


def get_evaporation(qr, tau_e, q_star, qv, qvs0, z):
    ev = qr / (tau_e * q_star) * np.max([approxfqv(z, qvs0) - qv, 0])
    return ev


def get_terminalvelocity(vt0, qr, q_star):
    vt = vt0 * qr / q_star
    # vt = 0
    return vt


def get_aerosolvelocity(vtnd, vt0, qr, q_star):
    vtn = vtnd + np.min([qr / q_star, 1]) * np.max([get_terminalvelocity(vt0, qr, q_star) - vtnd, 0])
    return vtn


def auxcdev(cte, qn0, qn, gamma, qv, qvs0, qr, tau_e, q_star, z):
    cdev = get_condensation(cte, qn0, qn, gamma, qv, qvs0, z) - get_evaporation(qr, tau_e, q_star, qv, qvs0, z)
    return cdev


def one_sided(dt, dz, u, vpar1, vpar2, space, tau_w):
    m = np.shape(u)[0]
    aux = np.zeros((m, 5))
    dzt = dt / dz
    vt0 = vpar2[0]
    vtnd = vpar2[1]
    q_star = vpar2[2]
    vel = u[:, 0]
    tem = u[:, 1]
    qv = u[:, 2]
    qr = u[:, 3]
    qn = u[:, 4]

    cflvt = np.zeros(m)
    cflvn = np.zeros(m)

    for i in range(1, m - 1):

        cflvt[i] = vel[i] - get_terminalvelocity(vt0, qr[i], q_star)
        cflvn[i] = vel[i] - get_aerosolvelocity(vtnd, vt0, qr[i], q_star)

        if vel[i] < 0:
            aux[i, 0] = vel[i] + dt * (get_bouyancyforce(vpar1, space[i], tem[i], qv[i], qr[i]) - vel[i] / tau_w)
            aux[i, 1] = tem[i] - dzt * vel[i] * (tem[i + 1] - tem[i])
            aux[i, 2] = qv[i] - dzt * vel[i] * (qv[i + 1] - qv[i])
        else:
            aux[i, 0] = vel[i] + dt * (get_bouyancyforce(vpar1, space[i], tem[i], qv[i], qr[i]))
            aux[i, 1] = tem[i] - dzt * vel[i] * (tem[i] - tem[i - 1])
            aux[i, 2] = qv[i] - dzt * vel[i] * (qv[i] - qv[i - 1])
        if cflvt[i] < 0:
            aux[i, 3] = qr[i] - dzt * ((vel[i] - get_terminalvelocity(vt0, qr[i + 1], q_star)) * qr[i + 1]
                                       - (cflvt[i] * qr[i]))
        else:
            aux[i, 3] = qr[i] - dzt * ((vel[i] - get_terminalvelocity(vt0, qr[i], q_star)) * qr[i] -
                                       (vel[i] - get_terminalvelocity(vt0, qr[i - 1], q_star)) * qr[i - 1])
        if cflvn[i] < 0:
            aux[i, 4] = (qn[i] - dzt *
                         ((vel[i] - get_aerosolvelocity(vtnd, vt0, qr[i + 1], q_star)) * qn[i + 1] - (
                                 vel[i] - get_aerosolvelocity(vtnd, vt0, qr[i], q_star)) * qn[i]))
        else:
            aux[i, 4] = (qn[i] - dzt *
                         ((vel[i] - get_aerosolvelocity(vtnd, vt0, qr[i], q_star)) * qn[i] - (
                                 vel[i - 1] - get_aerosolvelocity(vtnd, vt0, qr[i - 1], q_star)) * qn[i - 1]))
        aux[0, :] = aux[1, :]  # Antes de cambiarlo, funciona.
        aux[-1, :] = aux[-2, :]
        # print(get_terminalvelocity(vt0, aux[i, 3], q_star))
    pt = np.max(np.abs(cflvt))
    pn = np.max(np.abs(cflvn))
    p0 = np.max(np.abs(aux[:, 0]))
    p = np.max([pn, pt, p0])
    return aux, p

def g(qn0,vpar1,vpar2,vpar3,Z):
    tau_0 = vpar3[0]
    gamma = vpar3[1]
    q_star = vpar2[2]

    # y = auxcdev(tau_0,qn0,qn, gamma,qv=,qr=,tau_e=,q_star,Z)
    y = tau_0
    return y

def f(t):
    y = np.sin(t)
    return y


def one_sided(dt, dz, u, vpar1, vpar2, space, tau_w):
    m = np.shape(u)[0]
    aux = np.zeros((m, 5))
    dzt = dt / dz
    vt0 = vpar2[0]
    vtnd = vpar2[1]
    q_star = vpar2[2]
    vel = u[:, 0]
    tem = u[:, 1]
    qv = u[:, 2]
    qr = u[:, 3]
    qn = u[:, 4]

    cflvt = np.zeros(m)
    cflvn = np.zeros(m)

    for i in range(1, m - 1):

        cflvt[i] = vel[i] - get_terminalvelocity(vt0, qr[i], q_star)
        cflvn[i] = vel[i] - get_aerosolvelocity(vtnd, vt0, qr[i], q_star)

        if vel[i] < 0:
            aux[i, 0] = vel[i] + dt * (get_bouyancyforce(vpar1, space[i], tem[i], qv[i], qr[i]) - vel[i] / tau_w)
            aux[i, 1] = tem[i] - dzt * vel[i] * (tem[i + 1] - tem[i]) + dt * f(vel[i])
            aux[i, 2] = qv[i] - dzt * vel[i] * (qv[i + 1] - qv[i]) + dt * f(vel[i])
        else:
            aux[i, 0] = vel[i] + dt * (get_bouyancyforce(vpar1, space[i], tem[i], qv[i], qr[i]))
            aux[i, 1] = tem[i] - dzt * vel[i] * (tem[i] - tem[i - 1]) + dt * f(vel[i])
            aux[i, 2] = qv[i] - dzt * vel[i] * (qv[i] - qv[i - 1]) + dt * f(vel[i])
        if cflvt[i] < 0:
            aux[i, 3] = qr[i] - dzt * ((vel[i] - get_terminalvelocity(vt0, qr[i + 1], q_star)) * qr[i + 1]
                                       - (cflvt[i] * qr[i])) + dt * f(vel[i])
        else:
            aux[i, 3] = qr[i] - dzt * ((vel[i] - get_terminalvelocity(vt0, qr[i], q_star)) * qr[i] -
                                       (vel[i] - get_terminalvelocity(vt0, qr[i - 1], q_star)) * qr[i - 1]) + dt * f(vel[i])
        if cflvn[i] < 0:
            aux[i, 4] = (qn[i] - dzt *
                         ((vel[i] - get_aerosolvelocity(vtnd, vt0, qr[i + 1], q_star)) * qn[i + 1] - (
                                 vel[i] - get_aerosolvelocity(vtnd, vt0, qr[i], q_star)) * qn[i])) + dt * f(vel[i])
        else:
            aux[i, 4] = (qn[i] - dzt *
                         ((vel[i] - get_aerosolvelocity(vtnd, vt0, qr[i], q_star)) * qn[i] - (
                                 vel[i - 1] - get_aerosolvelocity(vtnd, vt0, qr[i - 1], q_star)) * qn[i - 1])) + dt * f(vel[i])
        aux[0, :] = aux[1, :]  # Antes de cambiarlo, funciona.
        aux[-1, :] = aux[-2, :]
        # print(get_terminalvelocity(vt0, aux[i, 3], q_star))
    pt = np.max(np.abs(cflvt))
    pn = np.max(np.abs(cflvn))
    p0 = np.max(np.abs(aux[:, 0]))
    p = np.max([pn, pt, p0])
    return aux, p




def resol_test(t_0, t_f, cfl, dt, dz, workspace, q_parameters, tick, basic_parameters, space, tau_w):
    tick_time = 0
    while t_0 < t_f:
        print(t_0)
        q = one_sided(dt, dz, workspace, basic_parameters, q_parameters, space, tau_w)
        workspace = q[0]
        dt = cfl * dz / np.abs(q[1])
        print("cfl = ", dt / dz * q[1])
        tick_time += 1
        if tick_time == tick:
            # print(workspace)
            break
        t_0 += dt
        # print(q[1])
    return workspace

