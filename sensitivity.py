"""
Select the least known parameter from literature (layer thickness or composition)
and study the sensitivity of that parameter with respect to Mass and Moment of Inertia
"""

import numpy as np
import matplotlib.pyplot as plt
import model_1 as m1

change = 20 #[%]
n = 10000

t_layers_l = np.load("t_layers_Low moi.npy")
t_layers_h = np.load("t_layers_High moi.npy")
t_layers_m = np.load("t_layers_Mean moi.npy")
rho_layers_l = np.load("rho_layers_Low moi.npy")
rho_layers_h = np.load("rho_layers_High moi.npy")
rho_layers_m = np.load("rho_layers_Mean moi.npy")

labels = np.array(["l", "h", "m"])

m_l = []
m_h = []
m_m = []
moi_l = []
moi_h = []
moi_m = []
t_l = []
t_h = []
t_m = []

for i in range(3):
    rho_layers = np.zeros(4)
    t_layers = np.zeros(4)
    t_list = []
    exec(f"rho_layers = rho_layers_{labels[i]}")
    exec(f"t_layers = t_layers_{labels[i]}")

    for t in np.arange(t_layers[-1] * (1 - change*1e-2), t_layers[-1] * (1 + change*1e-2), 1e3):
        t_layers_change = np.copy(t_layers)
        t_layers_change[-1] = t
        _, m, moi = m1.mass(t_layers_change, rho_layers, n)
        exec(f"moi_{labels[i]}.append(moi[-1])")
        exec(f"m_{labels[i]}.append(m[-1])")
        exec(f"t_{labels[i]}.append(t)")

    exec(f"plt.plot(t_{labels[i]}, m_{labels[i]})")
    plt.show()
    exec(f"plt.plot(t_{labels[i]}, moi_{labels[i]})")
    plt.show()