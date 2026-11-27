import numpy as np
import matplotlib.pyplot as plt
import model_1 as m1
import scipy as sp

def opt_1(m_ref, moi_ref, R_ref, t_layers, rho_layers, n):

    x_0 = np.concatenate([t_layers, rho_layers])
    x_norms = 10**(np.floor(np.log10(x_0)))
    x_0_scaled = x_0 / x_norms
    x_l = x_0_scaled * 0.9
    x_u = x_0_scaled * 1.1
    x_bounds = sp.optimize.Bounds(x_l, x_u)
    def obj(x_0_scaled):
        return m1.opt_obj(x_0_scaled, x_norms, n, R_ref, m_ref, moi_ref)
    sol = sp.optimize.minimize(obj, x_0_scaled, method='trust-constr', bounds=x_bounds)

    return sol.x * x_norms

if __name__ == "__main__":
    # Reference values
    m_ref = 4.8685e24  # [kg]
    moi_ref_l = 0.313
    moi_ref_h = 0.361
    moi_ref_m = (moi_ref_l + moi_ref_h) / 2
    R_ref = 6052e3  # [m]

    t_init = np.array([3228e3, 2091e3, 703e3, 30e3])       # [m]
    rho_init = np.array([10401.0, 4898.0, 3691.0, 2900.0])  # [kg/m^3]
    n = 10000

    mois_ref = np.array([moi_ref_l, moi_ref_h, moi_ref_m])
    labels = ["Low moi", "High moi", "Mean moi"]

    for i in range(mois_ref.size):
        x_opt = opt_1(m_ref, mois_ref[i], R_ref, t_init, rho_init, n)
        t_opt = x_opt[:4]
        rho_opt = x_opt[4:]
        r_M, M, moi = m1.mass(t_opt, rho_opt, n)
        print(f"=== {labels[i]} (target = {mois_ref[i]:.3f}) ===")
        print(f"t_layers [km]:      {np.round(t_opt * 1e-3, 2)} (R_out = {np.sum(t_opt)/1e3:.2f} km)")
        print(f"rho_layers [kg/m3]: {np.round(rho_opt, 1)}")
        print(
            f"Rel errors -> M: {(M[-1]-m_ref)/m_ref:.2e}, "
            f"MoI: {(moi[-1]-mois_ref[i])/mois_ref[i]:.2e}, "
            f"R: {(r_M[-1]-R_ref)/R_ref:.2e}\n"
        )
        np.save(f"t_layers_{labels[i]}", t_opt)
        np.save(f"rho_layers_{labels[i]}", rho_opt)