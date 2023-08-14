def initial_condition_w(x, t):
    if x < 0.5:
        y = 2
    else:
        y = 1
    return y


def initial_condition_theta(x, t):
    if x < 0.3:
        y = 1
    else:
        y = 0
    return y


def initial_condition_qv(x, t):
    if x < 0.5:
        y = 1
    else:
        y = 2
    return y


def initial_condition_qr(x, t):
    if x < 0.5:
        y = 0.5
    else:
        y = 1.5
    return y


def initial_condition_qn(x, t):
    if x < 0.5:
        y = 2
    else:
        y = 3
    return y
