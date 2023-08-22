import numpy as np

# Escalas
length_scale = 10  # Km
temperature_scale = 3  # K
time_scale = 15  # min
ratio_scale = 1000  # g/Kg Qs
velocity_scale = (length_scale * 1000) / (time_scale * 60)

# Parametros
B = 3  # K / km
B = B * length_scale / temperature_scale

heat_latence = 2.5e6
cp = 1e3
cp = cp * temperature_scale
lcp = heat_latence / cp

q_star = 1000000
q_star = q_star / ratio_scale

epsilon = 0.6

tau_w = 7.5  # min
tau_w = tau_w / time_scale

g = 9.81  # m/s^2
g = g / velocity_scale * 60 * time_scale

b_w = 0.491  # m/s sqrt(s)
b_w = b_w / (length_scale * 1000) * (time_scale * 60) * np.sqrt(time_scale * 60)

tau0 = 0.1
gamma = 0.3
taue = 0.2
c = 0.2
c0 = 0.4


vt0 = 1
vtnd = 1
