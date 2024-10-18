import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


if __name__ == "__main__":
    fig = plt.figure(figsize=(8, 5.5), layout="constrained")

    fs = "xx-large"
    subfigs = fig.subfigures(2, 1, height_ratios=[1, 1.3])

    ax_snaps = subfigs[0].subplots(1, 4)

    snaps = ["fig/snap/L800_e0.3_r1_rB0.1_n4.png",
             "fig/snap/L800_e0.3_r1_rB0.1_n2.png",
             "fig/snap/L800_e0.3_r1_rB0.1_n1.png",
             "fig/snap/L800_e0.3_r1.6_rB0.1_n1.png"
             ]
    
    labels = [
        r"(a) $\bar{\rho}_A=0.9$",
        r"(b) $\bar{\rho}_A=0.9$",
        r"(c) $\bar{\rho}_A=0.9$",
        r"(d) $\bar{\rho}_A=1.5$"
    ]
    for i, ax in enumerate(ax_snaps):
        im = mpimg.imread(snaps[i])
        ax.imshow(im)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        ax.set_title(labels[i], fontsize=fs)


    (ax1, ax2, ax3) = subfigs[1].subplots(1, 3, width_ratios=[1, 1, 0.9])

    profiles = ["data/profile/L800_800_d0.1000_e0.300_r1_s2004.npz",
                "data/profile/L800_800_d0.1000_e0.300_r1_s2002.npz",
                "data/profile/L800_800_d0.1000_e0.300_r1_s2011.npz"]
    labels = ["(a)", "(b)", "(c)"]
    colors = ["tab:red", "tab:green", "tab:blue"]
    for i, profile in enumerate(profiles):
        with np.load(profile, "r") as data:
            x, rhoA, mA = data["x"], data["rhoA"], data["mA"]
        ax1.plot(x, rhoA, label=labels[i], c=colors[i])
    ax1.set_xlim(0, 800)
    ax1.set_ylim(ymax=4.2)
    ax1.set_xlabel(r"$x$", fontsize=fs)
    ax1.set_title(r"(e) $\langle \rho_A(\mathbf{r}, t)\rangle_y$", fontsize=fs)
    ax1.legend(frameon=False, ncols=3, handlelength=1, fontsize="large", columnspacing=1)


    rho0_arr = np.array([1.0, 1.2, 1.4, 1.6, 1.8])
    m_arr = np.zeros_like(rho0_arr)
    for i, rho0 in enumerate(rho0_arr):
        filename = f"data/time_ave_profile/L800_800_rB0.1_e0.3_r{rho0:.2f}.npz"
        with np.load(filename, "r") as data:
            x, rhoA, mA = data["x"], data["rhoA"], data["mA"]
        line, = ax2.plot(x, np.roll(rhoA, 40), label=r"$%g$" % (rho0-0.1))
        m_arr[i] = np.mean(mA)
        ax3.plot(rho0_arr[i] - 0.1, m_arr[i], "o", c=line.get_c(), fillstyle="none")
    
    ax2.set_title(r"(f) $\langle \rho_A(\mathbf{r}, t)\rangle_{y, t}$", fontsize=fs)
    ax2.set_xlim(0, 800)
    ax2.set_xlabel(r"$x$", fontsize=fs)
    ax3.set_xlabel(r"$\bar{\rho}_A$", fontsize=fs)

    ax3.set_title(r"(g) $\langle \mathbf{w}_A(\mathbf{r},t)\cdot \hat{\mathbf{e}}_x\rangle_{\mathbf{r}, t}$", fontsize=fs)
    ax2.legend(title=r"$\bar{\rho}_A=$", ncols=1, frameon=False, title_fontsize="large")
    ax2.set_ylim(ymax=4.2)
    # p = np.polyfit(rho0_arr, m_arr, deg=1)
    # x = np.linspace(0, 1.8)
    # y = p[0] * x + p[1]
    # ax2.axhline(-p[1]/p[0])
    # ax3.plot(x, y)
    # ax3.set_xlim(0)
    # ax3.set_ylim(0)
    # plt.show()
    plt.savefig("fig/FIG3.pdf", dpi=200)
    plt.close()
