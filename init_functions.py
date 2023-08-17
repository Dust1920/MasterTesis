import parameters as p


def initial_function(a, b, s, x, t):
    if x < s:
        y = a
    else:
        y = b
    return y


# initial_parameters_functions
s_velocity = 1
s_velocity = s_velocity / p.velocity_scale

s_theta = 3
s_theta = s_theta / p.temperature_scale

s_qv = 28
s_qv = s_qv / p.ratio_scale
s_qr = 22
s_qr = s_qr / p.ratio_scale
s_qn = 26
s_qn = s_qn / p.ratio_scale
