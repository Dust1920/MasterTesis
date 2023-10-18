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


def one_sided_v5(dt, dz, u, vpar2, v0):
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
    for i in range(1, m - 1):
        if v0 < 0:
            aux[i, 0] = v0
            aux[i, 1] = tem[i] - dzt * v0 * (tem[i + 1] - tem[i])
            aux[i, 2] = qv[i] - dzt * v0 * (qv[i + 1] - qv[i])
            aux[i, 3] = (qr[i] - dzt *
                         ((v0 - get_terminalvelocity(vt0, qr[i + 1], q_star)) * qr[i + 1] - (
                                 v0 - get_terminalvelocity(vt0, qr[i], q_star)) * qr[i]))
            aux[i, 4] = (qn[i] - dzt *
                         ((v0 - get_aerosolvelocity(vtnd, vt0, qr[i + 1], q_star)) * qn[i + 1] - (
                                 v0 - get_aerosolvelocity(vtnd, vt0, qr[i], q_star)) * qn[i]))
        else:
            aux[i, 0] = v0
            aux[i, 1] = tem[i] - dzt * v0 * (tem[i] - tem[i - 1])
            aux[i, 2] = qv[i] - dzt * v0 * (qv[i] - qv[i - 1])
            aux[i, 3] = (qr[i] - dzt *
                         ((v0 - get_terminalvelocity(vt0, qr[i], q_star)) * qr[i] - (
                                 v0 - get_terminalvelocity(vt0, qr[i - 1], q_star)) * qr[i - 1]))
            aux[i, 4] = (qn[i] - dzt *
                         ((v0 - get_aerosolvelocity(vtnd, vt0, qr[i], q_star)) * qn[i] - (
                                 v0 - get_aerosolvelocity(vtnd, vt0, qr[i - 1], q_star)) * qn[i - 1]))
        aux[0, :] = aux[1, :]  # Antes de cambiarlo, funciona.
        aux[-1, :] = aux[-2, :]
    return aux


def one_sided_v6(dt, dz, u, vpar2, w):
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

        cflvt[i] = w - get_terminalvelocity(vt0, qr[i], q_star)
        cflvn[i] = w - get_aerosolvelocity(vtnd, vt0, qr[i], q_star)

        if w < 0:
            aux[i, 0] = w
            aux[i, 1] = tem[i] - dzt * w * (tem[i + 1] - tem[i])
            aux[i, 2] = qv[i] - dzt * w * (qv[i + 1] - qv[i])
        else:
            aux[i, 0] = w
            aux[i, 1] = tem[i] - dzt * w * (tem[i] - tem[i - 1])
            aux[i, 2] = qv[i] - dzt * w * (qv[i] - qv[i - 1])
        if cflvt[i] < 0:
            aux[i, 3] = qr[i] - dzt * ((w - get_terminalvelocity(vt0, qr[i + 1], q_star)) * qr[i + 1]
                                       - (cflvt[i] * qr[i]))
        else:
            aux[i, 3] = qr[i] - dzt * ((w - get_terminalvelocity(vt0, qr[i], q_star)) * qr[i] -
                                       (w - get_terminalvelocity(vt0, qr[i - 1], q_star)) * qr[i - 1])
        if cflvn[i] < 0:
            aux[i, 4] = (qn[i] - dzt *
                         ((w - get_aerosolvelocity(vtnd, vt0, qr[i + 1], q_star)) * qn[i + 1] - (
                                 w - get_aerosolvelocity(vtnd, vt0, qr[i], q_star)) * qn[i]))
        else:
            aux[i, 4] = (qn[i] - dzt *
                         ((w - get_aerosolvelocity(vtnd, vt0, qr[i], q_star)) * qn[i] - (
                                 w - get_aerosolvelocity(vtnd, vt0, qr[i - 1], q_star)) * qn[i - 1]))
        aux[0, :] = aux[1, :]  # Antes de cambiarlo, funciona.
        aux[-1, :] = aux[-2, :]
        # print(get_terminalvelocity(vt0, aux[i, 3], q_star))
    pt = np.max(np.abs(cflvt))
    pn = np.max(np.abs(cflvn))
    p0 = np.max(np.abs(aux[:, 0]))
    p = np.max([pn, pt, p0])
    return aux, p