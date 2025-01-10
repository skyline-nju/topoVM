import numpy as np
import matplotlib.pyplot as plt
from kinetic_linear_instab import lin_diagram_w_color
from hydro_linear_instab import find_eat_thresh, find_eta_para

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"



if __name__ == "__main__":
    fig = plt.figure(figsize=(8, 6.8), layout="constrained")
    subfigs = fig.subfigures(1, 2, width_ratios=[0.9, 0.10])
    axes = subfigs[0].subplots(2, 2, sharex=True, sharey=True)

    fins = ["data/hydros/r1_q0.0010_n400_nt60.npz",
            "data/kinetic/r01_q0.0010_K2_n200_nt60.npz",
            "data/kinetic/r01_q0.0010_K20_n200_nt60.npz",
            "data/kinetic/r01_q0.0010_K100_n200_nt60.npz"
    ]

    x1, y1 = find_eat_thresh()
    x2, y2 = find_eta_para()
    for i, ax in enumerate(axes.flat):
        im = lin_diagram_w_color(fins[i], ax, True)
        ax.set_xlim(0, 0.5)
        ax.set_ylim(0, 0.5)
        ax.plot(x1, y1, "-", c="w", lw=2, label=r"$\eta_t$")
        ax.plot(x2, y2, "--", c="tab:cyan", lw=2, label=r"$\eta_\parallel$")
        # if i == 0:
    fs = "xx-large"
    axes[0, 0].legend(fontsize=22, frameon=True, loc=(0.54, 0.65), borderpad=0.2)


    
    axes[0, 0].set_ylabel(r"$\eta$", fontsize=fs)
    axes[1, 0].set_ylabel(r"$\eta$", fontsize=fs)
    axes[1, 0].set_xlabel(r"$\bar{\rho}_B$", fontsize=fs)
    axes[1, 1].set_xlabel(r"$\bar{\rho}_B$", fontsize=fs)
    cax = subfigs[1].add_axes([0.01, 0.2, 0.35, 0.6])

    cb = subfigs[1].colorbar(im, cax=cax)
    cax.set_yticks([0, np.pi/4, np.pi/2])
    cax.set_yticklabels(["$0$", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$"])
    cax.tick_params(axis='y', which='major', labelsize="xx-large")
    axes[0, 0].set_title(r"(a) Hydrodynamic level", fontsize=fs)
    axes[0, 1].set_title(r"(b) Kinetic level: $K=2$", fontsize=fs)
    axes[1, 0].set_title(r"(c) Kinetic level: $K=20$", fontsize=fs)
    axes[1, 1].set_title(r"(d) Kinetic level: $K=100$", fontsize=fs)

    cb.set_label("Most unstable direction", fontsize=20)

    plt.show()
    # plt.savefig("fig/FIG4.pdf", dpi=600)
    plt.close()