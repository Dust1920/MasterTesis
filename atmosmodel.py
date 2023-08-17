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
    b = g / theta_0 * (theta - theta_hat + epsilon * theta_0 * (qv - qv_hat) - theta_0 * qr)
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
    return vt


def get_aerosolvelocity(vtnd, vt0, qr, q_star):
    vtn = vtnd + np.min([qr / q_star, 1]) * np.max([get_terminalvelocity(vt0, qr, q_star), 0])
    return vtn


def auxcdev(cte, qn0, qn, gamma, qv, qvs0, qr, tau_e, q_star, z):
    cdev = get_condensation(cte, qn0, qn, gamma, qv, qvs0, z) - get_evaporation(qr, tau_e, q_star, qv, qvs0, z)
    return cdev


def one_sided_model(dspace, dtime, u, z_vector, vpar, tau_w, vpar2):
    m = np.shape(u)[0]
    n = np.shape(u)[1]
    aux = np.zeros((m, n))
    dzt = dtime / dspace
    # First Column Velocity 0
    # Second Column Temperature 1
    # Third Column QV 2
    # Fourth Column QR 3
    # Fiveth Column QN 4
    vt0 = vpar2[0]
    vtnd = vpar2[1]
    q_star = vpar2[2]
    # w = -1
    # vt = 0
    # theta 0/3 1 3 > 0
    for i in np.arange(1, m - 1):
        aux[i, 0] = -1  # u[i, 0] + dtime * (
        # get_bouyancyforce(vpar, z_vector[i], u[i, 1], u[i, 2], u[i, 3]) - u[i, 0] / tau_w)
        aux[i, 1] = u[i, 1] - dzt * (u[i, 0] * u[i, 1] - u[i - 1, 0] * u[i - 1, 1])
        aux[i, 2] = u[i, 2] - dzt * (u[i, 0] * u[i, 2] - u[i - 1, 0] * u[i - 1, 2])
        aux[i, 3] = u[i, 3] - dzt * (
                (u[i, 0] - 0 * get_terminalvelocity(vt0, aux[i, 3], q_star)) * u[i, 3] -
                (u[i - 1, 0] - 0 * get_terminalvelocity(vt0, aux[i - 1, 3], q_star)) * u[i - 1, 3])
        aux[i, 4] = u[i, 4] - dzt * (
                (u[i, 0] - get_aerosolvelocity(vtnd, vt0, aux[i, 3], q_star)) * u[i, 4] -
                (u[i - 1, 0] - get_aerosolvelocity(vtnd, vt0, aux[i - 1, 3], q_star)) * u[i - 1, 4])
    return aux


def theta_0(z, b, t0):
    y = (t0 + b * z) * 0
    return y


def one_sided_model_2(dspace, dtime, u, vpar2, b, t0, qv0):
    m = np.shape(u)[0]
    n = np.shape(u)[1]
    aux = np.zeros((m, n))
    dzt = dtime / dspace
    # First Column Velocity 0
    # Second Column Temperature 1
    # Third Column QV 2
    # Fourth Column QR 3
    # Fiveth Column QN 4
    vt0 = vpar2[0]
    vtnd = vpar2[1]
    q_star = vpar2[2]
    for i in np.arange(1, m - 1):
        aux[i, 0] = -1
        aux[i, 1] = (u[i, 1] - dzt *
                     (u[i, 0] * (u[i, 1] + theta_0(i * dspace, b, t0)) - u[i - 1, 0] * (u[i - 1, 1] + theta_0(i * dspace, b, t0))))
        aux[i, 2] = (u[i, 2] - dzt * (u[i, 0] * (u[i, 2] + 0 * approxfqv(i * dspace, qv0)) - u[i - 1, 0] *
                                      (u[i - 1, 2] + 0 * approxfqv(i * dspace, qv0))))
        aux[i, 3] = u[i, 3] - dzt * (u[i, 0] * u[i, 3] - u[i - 1, 0] * u[i - 1, 3])
        aux[i, 4] = u[i, 4] - dzt * (u[i, 0] * u[i, 4] - u[i - 1, 0] * u[i - 1, 4])
    return aux