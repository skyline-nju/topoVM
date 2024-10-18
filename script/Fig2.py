import numpy as np
import matplotlib.pyplot as plt
from order_para import read_eta_phi_G_Grho_Gphi

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


if __name__ == "__main__":
    fig, axes = plt.subplots(3, 2, figsize=(8, 6), constrained_layout=True, sharey="row")

    dx = 40
    rhoB_arr = [0.03, 0.1]
    for j, rhoB in enumerate(rhoB_arr):
        if rhoB == 0.03:
            L_arr = [200, 400, 800]
        else:
            L_arr = [200, 400]
        for L in L_arr:
            eta, phi, G, G_rho, G_phi = read_eta_phi_G_Grho_Gphi(L, rhoB, dx=dx)
            axes[0, j].plot(eta, phi, "-o", fillstyle="none")
            axes[1, j].plot(eta, G_rho, '-o', fillstyle="none")
    

    axes[1, 0].set_ylim(ymin=-1, ymax=0.7)
    axes[1, 0].axhline(2/3, linestyle="dashed", color="tab:grey")
    axes[1, 1].axhline(2/3, linestyle="dashed", color="tab:grey")
    # plt.show()
    plt.savefig("fig/FIG2.pdf")
    plt.close()
