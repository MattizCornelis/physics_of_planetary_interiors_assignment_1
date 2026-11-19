"""
Develop a multiple layered homogenous density structure that uses average densities
from the composition you. Calculate your own internal pressure and gravity profiles.
Check if your model fit with the observed: Mass, Radius and Moment of Inertia. Selected
the least known parameter from literature (layer thickness or composition) and study the
sensitivity of that parameter with respect to Mass and Moment of Inertia. Make 3 models
that fit the observations (min/mean/max)
"""

import numpy as np
import matplotlib.pyplot as plt

def dM(r, rho):
    return 4 * np.pi * rho * r ** 2

def Euler(M_0, r_0, dr, rho, dM):
    r = r_0 + dr
    M = M_0 + dr * dM(r_0, rho)
    return r, M

def rho(rho_layers, t_layers, r):
    idx = np.where(r <= r_layers(t_layers))[0][0]
    return rho_layers[idx]

def r_layers(t_layers):
    r_layers = np.zeros_like(t_layers)
    for i in range(t_layers.shape[0]):
        r_layers[i] = t_layers[i] + r_layers[i-1]
    return r_layers

def mass(t_layers, rho_layers, n):
    r_0 = 0
    M_0 = 0
    R = np.sum(t_layers)
    dr = R / n
    mass_log = [M_0]
    r_log = []
    for r in np.arange(r_0, R, dr):
        r, M = Euler(mass_log[-1], r, dr, rho(rho_layers, t_layers, r), dM)
        r_log.append(r)
        mass_log.append(M)

    mass_log.pop(0)

    return r_log, mass_log

def gravity(M, r, n):
    G = 6.674e-11
    r_g = r.copy()
    r_g.pop(0)
    r_g = np.array(r_g)
    M_g = M.copy()
    M_g.pop(0)
    M_g = np.array(M_g)
    print(r[0:3], r_g[0:2], M[0:3], M_g[0:2])
    return G * M_g / (r_g)**2

if __name__ == '__main__':
    t_layers = np.array([1, 1, 1])
    rho_layers = np.array([3, 2, 1])
    n = 1e3
    r, M = mass(t_layers, rho_layers, n)
    g = gravity(r, M, n)
    plt.plot(r, M)
    plt.show()
    plt.plot(r[1:], g)
    plt.show()