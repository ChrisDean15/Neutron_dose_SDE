import math
import numpy as np
import random


def log_a(k, m, theta):
    ret = (
        math.log(theta + 2 * k - 1)
        + math.lgamma(theta + m + k - 1)
        - math.lgamma(theta + m)
        - math.lgamma(m + 1)
        - math.lgamma(k - m + 1)
    )
    return ret


def b(k, m, t, theta):
    ret = 1
    if k > 0:
        ret = math.exp(log_a(k, m, theta) - k * (k + theta - 1) * t / 2)
    return ret


def c(m, t, theta):
    i = 0
    b_curr = b(m, m, t, theta)
    b_next = b(m + 1, m, t, theta)
    while b_next >= b_curr:
        i += 1
        b_curr = b_next
        b_next = b(i + m + 1, m, t, theta)
    return i


def number_of_blocks(t):
    m = 0
    theta = 1
    if t < 0.07:
        mu = 2 / t
        sigma = math.sqrt(2 / (3 * t))
        m = round(random.gauss(mu, sigma))
    else:
        k = [0]
        proceed = True
        u = random.random()
        smin = 0
        smax = 0
        increment = 0
        while proceed:
            k[m] = math.ceil(c(m, t, theta) / 2)
            for i in range(k[m]):
                increment = b(m + 2 * i, m, t, theta) - b(m + 2 * i + 1, m, t, theta)
                smin += increment
                smax += increment
            increment = b(m + 2 * k[m], m, t, theta)
            smin += increment - b(m + 2 * k[m] + 1, m, t, theta)
            smax += increment
            while smin < u and u < smax:
                for i in range(m + 1):
                    k[i] += 1
                    increment = b(i + 2 * k[i], i, t, theta)
                    smax = smin + increment
                    smin += increment - b(i + 2 * k[i] + 1, i, t, theta)
            if smin > u:
                proceed = False
            else:
                k.append(0)
                m += 1
    return m


def wright_fisher(r):
    if r > 1e-9:
        m = number_of_blocks(r)
        y = np.random.beta(1, 1 + m)
    else:
        y = r / 2
        y = fabs(random.gauss(0, sqrt(r * y * (1 - y))))
    return y
