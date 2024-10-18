import numpy as np
import matplotlib.pyplot as plt
from kinetic_linear_instab import plot_rhoA_eta_panel, plot_rhoB_eta_panel

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"



if __name__ == "__main__":
    fig, axes = plt.subplots(2, 2, figsize=(8, 6), constrained_layout=True, sharey="row", sharex="col")

    fs = "xx-large"


    plot_rhoB_eta_panel(axes[1, 0], 100, 200)
    plot_rhoA_eta_panel(axes[1, 1], 100, 200)
    axes[0, 0].set_ylabel(r"$\eta$", fontsize=fs)
    axes[1, 0].set_ylabel(r"$\eta$", fontsize=fs)
    axes[1, 0].set_xlabel(r"$\bar{\rho}_B$", fontsize=fs)
    axes[1, 1].set_xlabel(r"$\bar{\rho}_A$", fontsize=fs)
    # plt.show()
    plt.savefig("fig/FIG4.pdf")
    plt.close()