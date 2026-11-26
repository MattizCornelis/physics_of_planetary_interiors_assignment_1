import numpy as np
import matplotlib.pyplot as plt
import model_1 as m1

def opt_1(m_ref, moi_ref_l, moi_ref_h, moi_ref_m, R_ref, t_layers, rho_layers, n):
    r_M, M, moi = m1.mass(t_layers, rho_layers, n)
    g = m1.gravity(M, r_M)

    m_out = M[-1]
    moi_out = moi[-1]
    R_out = r_M[-1]

    m_err = (m_out - m_ref) / m_ref
    moi_err_l = (moi_out - moi_ref_l) / moi_ref_l
    moi_err_h = (moi_out - moi_ref_h) / moi_ref_h
    moi_err_m = (moi_out - moi_ref_m) / moi_ref_m
    R_err = (R_out - R_ref) / R_ref


if __name__ == "__main__":
    # Reference values

    m_ref = 4.8685e24  # [kg]
    moi_ref_l = 0.313
    moi_ref_h = 0.361
    moi_ref_m = (moi_ref_l + moi_ref_h) / 2
    R_ref = 6052e3  # [km]

    t_layers = np.array([3228e3, 2091e3, 703e3, 30e3])  # [m]
    rho_layers = np.array([10401, 4898, 3691, 2900])  # [kg/m^3]
    n = 1000