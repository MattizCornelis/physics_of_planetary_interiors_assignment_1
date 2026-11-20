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

def dp(rho, g):
    return - rho * g

def dmoi(r, rho):
    return 8 * np.pi / 3 * rho * r ** 4

def Euler_M(M_0, r_0, dr, rho, dM):
    r = r_0 + dr
    M = M_0 + dr * dM(r_0, rho)
    return r, M

def Euler_p(p_0, r_0, dr, rho, g, dp):
    r = r_0 - dr
    p = p_0 - dr * dp(rho, g)
    return r, p

def Euler_moi(moi_0, r_0, dr, rho, dmoi):
    r = r_0 + dr
    moi = moi_0 + dr * dmoi(r_0, rho)
    return r, moi

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
    moi_0 = 0
    R = np.sum(t_layers)
    dr = R / n
    mass_log = [M_0]
    moi_log = [moi_0]
    r_log = []
    for r in np.arange(r_0, R, dr):
        r, M = Euler_M(mass_log[-1], r, dr, rho(rho_layers, t_layers, r), dM)
        _, moi = Euler_moi(moi_log[-1], r, dr, rho(rho_layers, t_layers, r), dmoi)
        r_log.append(r)
        mass_log.append(M)
        moi_log.append(moi)

    mass_log.pop(0)
    moi_log.pop(0)

    return r_log, mass_log, moi_log

def gravity(M, r):
    G = 6.674e-11
    r_g = r.copy()
    r_g.pop(0)
    r_g = np.array(r_g)
    M_g = M.copy()
    M_g.pop(0)
    M_g = np.array(M_g)
    return G * M_g / (r_g)**2

def pressure(t_layers, rho_layers, n, g_lst, r_lst):
    R = np.sum(t_layers)
    r_0 = R
    p_0 = 0
    dr = R / n
    pressure_log = [p_0]
    r_log = [r_0]
    for r in np.arange(r_0, 0, -dr):
        idx_g = np.where(r <= r_lst)[0][0]
        g = g_lst[idx_g]
        r, p = Euler_p(pressure_log[-1], r, dr, rho(rho_layers, t_layers, r), g, dp)
        r_log.append(r)
        pressure_log.append(p)
    pressure_log.pop(-1)
    r_log.pop(-1)

    return r_log, pressure_log

if __name__ == '__main__':
    t_layers = np.array([1, 1, 1, 1])
    rho_layers = np.array([4, 3, 2, 1])
    n = 1e3
    r_M, M, moi = mass(t_layers, rho_layers, n)
    g = gravity(r_M, M, n)
    r_p, p = pressure(t_layers, rho_layers, n, g, r_M[1:])
    plt.plot(r_M, M)
    plt.show()
    plt.plot(r_M, moi)
    plt.show()
    plt.plot(r_M[1:], g)
    plt.loglog()
    plt.show()
    plt.plot(r_p, p)
    plt.loglog()
    plt.show()