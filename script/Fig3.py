import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
from order_para import read_order_para_series
from GNF import read_GNF
from band_profile import plot_band_profiles_varied_rhoB, plot_band_profiles_varied_eta
from add_line import add_line

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


if __name__ == "__main__":
    fig = plt.figure(figsize=(8, 5), layout="constrained")
    subfigs = fig.subfigures(1, 2, width_ratios=[2.3, 1])
    ax_snap = subfigs[1].subplots(2, 1)
    axes = subfigs[0].subplots(2, 2, sharex=True, sharey="row")


    snaps = ["fig/snap/L400_e0.1_rB0.03.png",
            "fig/snap/L400_e0.1_rB0.30.png"
            # "fig/snap/L400_e0.3_rB0.30.png"
            ]

    for i, ax in enumerate(ax_snap):
            im = mpimg.imread(snaps[i])
            ax.imshow(im)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")
            # ax.set_title(labels[i], fontsize=fs)
    
    rhoB_arr = [0.03, 0.1, 0.2, 0.3]
    plot_band_profiles_varied_rhoB(axes[:, 0], rhoB_arr, lw=1.5)
    axes[0, 0].set_yscale("log")

    plot_band_profiles_varied_eta(axes[:, 1], lw=1.5)

    axes[0, 0].set_xlim(0, 400)

    axes[1, 1].set_xlabel(r"$x$", fontsize="xx-large")
    axes[1, 0].set_xlabel(r"$x$", fontsize="xx-large")
    axes[0, 0].set_ylabel(r"$\langle \rho_A \rangle_{y, t} $", fontsize="xx-large")
    # axes[1, 0].set_ylabel(r"$\frac{\langle m_{x,A}\rangle_{y, t}}{\langle \rho_A \rangle_{y, t}}$", fontsize="xx-large")
    axes[1, 0].set_ylabel(r"$\langle m_{x,A}\rangle_{y, t} \slash \langle \rho_A \rangle_{y, t}$", fontsize="xx-large")

    ax_snap[0].set_title(r"$\eta=0.1$", fontsize="xx-large")
    axes[0, 0].set_title(r"$\eta=0.1$", fontsize="xx-large")
    axes[0, 1].set_title(r"$\bar{\rho}_B=0.3$", fontsize="xx-large")

    axes[0, 0].legend(frameon=False, handlelength=1, title=r"$\bar{\rho}_B=$", labelspacing=0.2, title_fontsize="x-large", fontsize="large")
    axes[0, 1].legend(frameon=False, handlelength=1, title=r"$\eta=$", labelspacing=0.2, title_fontsize="x-large", fontsize="large")

    ax_snap[0].text(0.01, 0.9, r"(c) $\bar{\rho}_B=0.03$", transform=ax_snap[0].transAxes, fontsize="xx-large")
    ax_snap[1].text(0.01, 0.9, r"(f) $\bar{\rho}_B=0.3$", transform=ax_snap[1].transAxes, fontsize="xx-large")

    axes[0, 0].text(0.01, 0.88, r"(a)", transform=axes[0, 0].transAxes, fontsize="xx-large")
    axes[1, 0].text(0.01, 0.88, r"(d)", transform=axes[1, 0].transAxes, fontsize="xx-large")
    axes[0, 1].text(0.01, 0.88, r"(b)", transform=axes[0, 1].transAxes, fontsize="xx-large")
    axes[1, 1].text(0.01, 0.88, r"(e)", transform=axes[1, 1].transAxes, fontsize="xx-large")


    ax_snap[0].spines['bottom'].set_color("tab:blue")
    ax_snap[0].spines['top'].set_color("tab:blue")
    ax_snap[0].spines['left'].set_color("tab:blue")
    ax_snap[0].spines['right'].set_color("tab:blue")

    ax_snap[1].spines['bottom'].set_color("tab:red")
    ax_snap[1].spines['top'].set_color("tab:red")
    ax_snap[1].spines['left'].set_color("tab:red")
    ax_snap[1].spines['right'].set_color("tab:red")

    ax_snap[0].arrow(0.1, 0.13, 0.15, 0, transform=ax_snap[0].transAxes, color="k", head_length=0.04, width=0.012)
    ax_snap[0].arrow(0.1, 0.13, 0, 0.15, transform=ax_snap[0].transAxes, color="k", head_length=0.04, width=0.012)
    ax_snap[0].text(0.2, 0.03, r"$x$", transform=ax_snap[0].transAxes, fontsize="xx-large")
    ax_snap[0].text(0.01, 0.3, r"$y$", transform=ax_snap[0].transAxes, fontsize="xx-large")

    dx = 0.25
    ax_cb = ax_snap[1].inset_axes([0, 0, dx, dx])
    ax_cb.set_xticklabels([])
    ax_cb.set_yticklabels([])
    ax_cb.set_xticks([])
    ax_cb.set_yticks([])
    im = mpimg.imread("fig/circle2.png")
    ax_cb.imshow(im)
    # ax_cb.set_title(r"$\theta_i$", fontsize=fs)
    ax_cb.set_title(r"$\theta_i$", fontsize="xx-large")

    plt.show()
#     plt.savefig("fig/FIG3.pdf", dpi=200)
    plt.close()