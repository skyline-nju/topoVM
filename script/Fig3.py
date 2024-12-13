import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
from order_para import read_order_para_series
from GNF import read_GNF
from band_profile import plot_band_profiles_varied_rhoB, plot_band_profiles_varied_eta
from add_line import add_line

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


def six_panels():
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


def seven_panels():
    fig = plt.figure(figsize=(8, 7.5), layout="constrained")
    subfigs = fig.subfigures(2, 1, height_ratios=[1, 1.3], hspace=0.)
    ax_snap = subfigs[0].subplots(1, 3)
    axes = subfigs[1].subplots(2, 2, sharex=True, sharey="row")


    snaps = ["fig/snap/L400_e0.1_rB0.03.png",
            "fig/snap/L400_e0.1_rB0.30.png",
            "fig/snap/L400_e0.3_rB0.30.png"
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
    plot_band_profiles_varied_rhoB(axes[:, 0], rhoB_arr, lw=2)
    axes[0, 0].set_yscale("log")

    plot_band_profiles_varied_eta(axes[:, 1], lw=2)

    axes[0, 0].set_xlim(0, 400)

    axes[1, 1].set_xlabel(r"$x$", fontsize="xx-large")
    axes[1, 0].set_xlabel(r"$x$", fontsize="xx-large")
    axes[0, 0].set_ylabel(r"$\langle \rho_A \rangle_{y, t} $", fontsize="xx-large")
    # axes[1, 0].set_ylabel(r"$\frac{\langle m_{x,A}\rangle_{y, t}}{\langle \rho_A \rangle_{y, t}}$", fontsize="xx-large")
    axes[1, 0].set_ylabel(r"$\langle m_{x,A}\rangle_{y, t} \slash \langle \rho_A \rangle_{y, t}$", fontsize="xx-large")

#     ax_snap[0].set_title(r"$\eta=0.1$", fontsize="xx-large")
    axes[0, 0].set_title(r"$\eta=0.1$", fontsize="xx-large")
    axes[0, 1].set_title(r"$\bar{\rho}_B=0.3$", fontsize="xx-large")

    axes[0, 0].legend(frameon=False, handlelength=1, title=r"$\bar{\rho}_B=$", labelspacing=0.2, title_fontsize="x-large", fontsize="x-large", loc=(0.7, 0.21))
    axes[0, 1].legend(frameon=False, handlelength=1, title=r"$\eta=$", labelspacing=0.2, title_fontsize="x-large", fontsize="x-large", loc=(0.7, 0.35))

    ax_snap[0].text(0.01, 0.9, r"(a)", transform=ax_snap[0].transAxes, fontsize="xx-large")
    ax_snap[1].text(0.01, 0.9, r"(b)", transform=ax_snap[1].transAxes, fontsize="xx-large")
    ax_snap[2].text(0.01, 0.9, r"(c)", transform=ax_snap[2].transAxes, fontsize="xx-large")

    axes[0, 0].text(0.01, 0.88, r"(d)", transform=axes[0, 0].transAxes, fontsize="xx-large")
    axes[1, 0].text(0.01, 0.88, r"(f)", transform=axes[1, 0].transAxes, fontsize="xx-large")
    axes[0, 1].text(0.01, 0.88, r"(e)", transform=axes[0, 1].transAxes, fontsize="xx-large")
    axes[1, 1].text(0.01, 0.88, r"(g)", transform=axes[1, 1].transAxes, fontsize="xx-large")


    ax_snap[0].spines['bottom'].set_color("tab:blue")
    ax_snap[0].spines['top'].set_color("tab:blue")
    ax_snap[0].spines['left'].set_color("tab:blue")
    ax_snap[0].spines['right'].set_color("tab:blue")

    ax_snap[1].spines['bottom'].set_color("tab:red")
    ax_snap[1].spines['top'].set_color("tab:red")
    ax_snap[1].spines['left'].set_color("tab:red")
    ax_snap[1].spines['right'].set_color("tab:red")

    ax_snap[2].spines['bottom'].set_color("tab:brown")
    ax_snap[2].spines['top'].set_color("tab:brown")
    ax_snap[2].spines['left'].set_color("tab:brown")
    ax_snap[2].spines['right'].set_color("tab:brown")

    ax_snap[0].set_title(r"$(\bar{\rho}_B,\eta)=(0.03, 0.1)$", fontsize="xx-large", color="tab:blue")
    ax_snap[1].set_title(r"$(\bar{\rho}_B,\eta)=(0.3, 0.1)$", fontsize="xx-large", color="tab:red")
    ax_snap[2].set_title(r"$(\bar{\rho}_B,\eta)=(0.3, 0.3)$", fontsize="xx-large", color="tab:brown")

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


if __name__ == "__main__":
    fig, axes = plt.subplots(2, 2, constrained_layout=True, sharex="col", figsize=(8, 6), width_ratios=[1, 1.5])
    ax_snap = axes[:, 0]
    ax_profile = axes[:, 1]


    snaps = ["fig/snap/L400_e0.1_rB0.03.png",
            "fig/snap/L400_e0.1_rB0.30.png",
            ]

    for i, ax in enumerate(ax_snap):
            im = mpimg.imread(snaps[i])
            ax.imshow(im)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")
    
    rhoB_arr = [0.03, 0.1, 0.2, 0.3]
    plot_band_profiles_varied_rhoB(ax_profile, rhoB_arr, lw=2)
    ax_profile[0].set_yscale("log")

    ax_profile[0].set_xlim(0, 400)

    ax_profile[0].set_ylabel(r"$\langle \rho_A \rangle_{y, t} $", fontsize="xx-large")
    ax_profile[1].set_ylabel(r"$\langle m_{x,A}\rangle_{y, t} \slash \langle \rho_A \rangle_{y, t}$", fontsize="xx-large")
    ax_profile[1].set_ylim(-0.15)
    ax_profile[0].legend(frameon=False, handlelength=1, title=r"$\bar{\rho}_B=$", labelspacing=0.2, title_fontsize="xx-large", fontsize="x-large", loc=(0.05, 0.4))

    ax_in = ax_profile[0].inset_axes([0.65, 0.4, 0.35, 0.6])
    plot_band_profiles_varied_rhoB(ax_in, rhoB_arr, lw=1.5)
    ax_in.set_xlim(0, 400)
    ax_in.set_ylim(0, 14.5)
    ax_in.text(0.85, 0.05, r"$x$", transform=ax_in.transAxes, fontsize="xx-large")
    ax_in.text(0.03, 0.6, r"$\langle \rho_A \rangle_{y,t}$", transform=ax_in.transAxes, fontsize="xx-large", rotation=90)


    bbox = dict(edgecolor="w", facecolor="w", boxstyle="Square, pad=0.03")
    ax_snap[0].text(0.01, 0.915, r"(a) $\bar{\rho}_B=0.03$", transform=ax_snap[0].transAxes, fontsize="xx-large", bbox=bbox)
    ax_snap[1].text(0.01, 0.92, r"(b) $\bar{\rho}_B=0.3$", transform=ax_snap[1].transAxes, fontsize="xx-large", bbox=bbox)
    ax_profile[0].text(0.01, 0.92, r"(c)", transform=ax_profile[0].transAxes, fontsize="xx-large")
    ax_profile[1].text(0.01, 0.92, r"(d)", transform=ax_profile[1].transAxes, fontsize="xx-large")
    ax_profile[1].text(0.9, 0.03, r"$x$", transform=ax_profile[1].transAxes, fontsize="xx-large")
    ax_profile[0].text(0.9, 0.03, r"$x$", transform=ax_profile[0].transAxes, fontsize="xx-large")

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

#     plt.show()
    plt.savefig("fig/FIG3.pdf", dpi=200)
    plt.close()
