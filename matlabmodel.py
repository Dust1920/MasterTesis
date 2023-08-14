import parameters as p
import numpy as np
import atmosmodel as at
import matplotlib.pyplot as plt


def bouyancy_matlab(z, par):
    b = par[0]
    theta_0 = par[1]
    qv0 = par[2]
    qvs0 = par[3]
    lcp = par[4]
    epsilon = par[5]
    theta_e = par[6]
    qt = par[7]
    g = par[8]

    b1 = theta_0 + z * b
    b2 = at.approxfqv(z, qv0)
    b3 = at.approxfqv(z, qvs0)
    b4 = b1 + lcp * b2
    b5 = b4 + theta_0 * (epsilon - lcp / theta_0) * b2
    b6 = np.min([b3, qt])
    b7 = np.max([qt - b3, 0])
    b8 = theta_e + theta_0 * (epsilon - lcp / theta_0) * b6 - theta_0 * b7
    b = g / theta_0 * (b8 - b5)
    return b


# Condiciones iniciales
zi = 0
zi = zi / p.length_scale
zfinal = 15
zfinal = zfinal / p.length_scale
zpts = 101
dz = (zfinal - zi) / (zpts - 1)
height = np.arange(zi, zfinal, dz)

z_initial = 0.5
z_initial = z_initial / p.length_scale

theta_0 = 300
theta_0 = theta_0 / p.temperature_scale

qv_0 = 16
qv_0 = qv_0 / p.ratio_scale

qvs_0 = 28
qvs_0 = qvs_0 / p.ratio_scale

theta_e = (theta_0 + p.B * z_initial) + p.lcp * at.approxfqv(z_initial, qv_0)
qt = at.approxfqv(z_initial, qv_0)

vpar = np.array([p.B, theta_0, qv_0, qvs_0, p.lcp, p.epsilon, theta_e, qt, p.g])

f = np.array([bouyancy_matlab(i, vpar) for i in height])

plt.plot(f * p.velocity_scale, height * p.length_scale)
plt.show()