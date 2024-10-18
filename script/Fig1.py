import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
from sim_PD import rhoA_vs_eta, rhoB_vs_eta

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


if __name__ == "__main__":
    fig = plt.figure(figsize=(8, 7.1), layout="constrained")

    fs = "xx-large"
    subfigs = fig.subfigures(2, 1, height_ratios=[4, 2.8])

    (ax1, ax2) = subfigs[0].subplots(1, 2, sharey=True)

    rhoB_vs_eta(ax1)
    rhoA_vs_eta(ax2)

    ax1.set_title(r"(a) $\rho_0=1$", fontsize=fs)
    ax2.set_title(r"(b) $\bar{\rho}_B=0.03$", fontsize=fs)

    for tick in ax1.get_yticklabels():
        tick.set_rotation(90)

    bbox = dict(edgecolor="w", facecolor="w", boxstyle="Square, pad=0.05")
    ax1.text(0.01, 0.94, r"$\eta$", transform=ax1.transAxes, fontsize=fs, bbox=bbox)
    ax2.text(0.01, 0.94, r"$\eta$", transform=ax2.transAxes, fontsize=fs, bbox=bbox)

    ax1.text(0.85, 0.03, r"$\bar{\rho}_B$", transform=ax1.transAxes, fontsize=fs, bbox=bbox)
    ax2.text(0.90, 0.03, r"$\bar{\rho}_A$", transform=ax2.transAxes, fontsize=fs, bbox=bbox)

    ax1.plot(0.03, 0.1, "v", c="tab:pink", ms=8)
    ax1.plot(0.03, 0.3, "p", c="tab:green", ms=8)
    ax1.plot(0.03, 0.42, "^", c="tab:red", ms=8)
    # ax1.set_ylabel(r"$\eta$", fontsize=fs)
    # ax1.set_xlabel(r"$\bar{\rho}_B$", fontsize=fs)
    # ax2.set_xlabel(r"$\bar{\rho}_A$", fontsize=fs)
    ax2.legend(fontsize="x-large", loc="lower center", framealpha=0.88, markerscale=1.5)
    ax_snaps = subfigs[1].subplots(1, 3)

    snaps = ["fig/snap/L2000_rB0.03_e0.10.png",
            "fig/snap/L2000_rB0.03_e0.30.png",
            "fig/snap/L2000_rB0.03_e0.42.png"]
    labels = [
        r"(c) $\eta=0.1$",
        r"(d) $\eta=0.3$",
        r"(e) $\eta=0.42$",
    ]

    mk = ["v", "p", "^"]
    mkc = ["tab:pink", "tab:green", "tab:red"]
    for i, ax in enumerate(ax_snaps):
        im = mpimg.imread(snaps[i])
        ax.imshow(im)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        ax.set_title(labels[i], fontsize=fs)
    
        dx = 0.1
        ax_in = ax.inset_axes([1-dx, 1-dx, dx, dx])
        ax_in.set_xticklabels([])
        ax_in.set_yticklabels([])
        ax_in.set_xticks([])
        ax_in.set_yticks([])
        ax_in.plot(0, 0, mk[i], c=mkc[i], ms=8)
    
    dx = 0.2
    ax_cb = ax_snaps[1].inset_axes([1-dx, 0, dx, dx])
    ax_cb.set_xticklabels([])
    ax_cb.set_yticklabels([])
    ax_cb.set_xticks([])
    ax_cb.set_yticks([])
    im = mpimg.imread("fig/circle2.png")
    ax_cb.imshow(im)
    ax_cb.set_title(r"$\theta_i$", fontsize=fs)
    
    # plt.show()
    plt.savefig("fig/FIG1.pdf", dpi=200)
    plt.close()
    


