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
    r_layers = np.cumsum(t_layers)
    if r >= r_layers[-1]:
        return rho_layers[-1]
    idx = np.where(r < r_layers)[0][0]
    return rho_layers[idx]

def mass(t_layers, rho_layers, n):
    R = np.sum(t_layers)
    dr = R / n

    M_curr = 0
    moi_curr = 0

    mass_log = []
    moi_log = []
    r_log = []

    for r_curr in np.arange(0, R, dr):
        current_rho = rho(rho_layers, t_layers, r_curr)

        r_next, M_next = Euler_M(M_curr, r_curr, dr, current_rho, dM)
        _, moi_next = Euler_moi(moi_curr, r_curr, dr, current_rho, dmoi)

        r_log.append(r_next)
        mass_log.append(M_next)
        moi_log.append(moi_next)

        M_curr = M_next
        moi_curr = moi_next

    for i in range(len(moi_log)):
        if r_log[i] != 0.0 and mass_log[i] != 0.0:
            moi_log[i] /= mass_log[i] * r_log[i]**2

    return r_log, mass_log, moi_log

def gravity(M, r):
    G = 6.674e-11
    return G * np.array(M) / (np.array(r)**2)

def pressure(t_layers, rho_layers, n, g_lst, r_lst):
    R = np.sum(t_layers)
    dr = R / n
    p_current = 0

    pressure_log = [p_current]
    r_log = [R]

    for r_curr in np.arange(R, 0, -dr):
        current_rho = rho(rho_layers, t_layers, r_curr)

        idx_g = np.argmin(np.abs(np.array(r_lst) - r_curr))
        g = g_lst[idx_g]

        r_next, p_next = Euler_p(p_current, r_curr, dr, current_rho, g, dp)

        r_log.append(r_next)
        pressure_log.append(p_next)
        p_current = p_next

    return r_log, pressure_log

if __name__ == '__main__':
    t_layers = np.array([3228e3, 2091e3, 703e3, 30e3]) # [m]
    rho_layers = np.array([10401, 4898, 3691, 2900]) # [kg/m^3]
    n = 1000

    r_M, M, moi = mass(t_layers, rho_layers, n)
    g = gravity(M, r_M)
    r_p, p = pressure(t_layers, rho_layers, n, g, r_M)

    for i in range(len(r_M)):
        r_M[i] = r_M[i] * 1e-3

    for i in range(len(r_p)):
        r_p[i] = r_p[i] * 1e-3
        p[i] = p[i] * 1e-9

    print("--- Model Output Values ---")
    print(f"Total Mass: {M[-1]:.4e} kg")
    print(f"Total Mass Moment of Inertia: {moi[-1]:.4}")
    print(f"Total Radius: {r_M[-1]} km")
    print(f"Total Surface Gravity: {g[-1]:.4} m/s^2")

    plt.plot(r_M, M)
    plt.xlabel('Radius [km]')
    plt.ylabel('Mass [kg]')
    plt.show()

    plt.plot(r_M, moi)
    plt.xlabel('Radius [km]')
    plt.ylabel('Moment of Inertia [-]')
    plt.show()

    plt.plot(r_M, g)
    plt.xlabel('Radius [km]')
    plt.ylabel(r'Surface Gravity [$\frac{m}{s^2}$]')
    plt.show()

    plt.plot(r_p, p)
    plt.xlabel('Radius [km]')
    plt.ylabel('Pressure [GPa]')
    plt.show()