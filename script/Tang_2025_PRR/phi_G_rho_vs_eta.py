import numpy as np
import matplotlib.pyplot as plt
from order_para import read_eta_phi_G_Grho_Gphi, get_order_para_varied_dx

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


def plot_phi_G_rho_vs_eta(axes, rhoB, dx=40):
    mk = {200:"-o", 400:"-s", 800: "->"}
    clist = {200: "tab:blue", 400: "tab:orange", 800: "tab:green"}

    if rhoB == 0.03:
        L_arr = [200, 400, 800]
    else:
        L_arr = [200, 400, 800]
    for L in L_arr:
        eta, phi, G, G_rho, G_phi = read_eta_phi_G_Grho_Gphi(L, rhoB, dx=dx)
        # if j == 0:
        #     label = None
        # else:
        #     label=r"$L=%g$" % L
        axes[0].plot(eta, phi, mk[L], fillstyle="none", ms=5, lw=1, c=clist[L], label=r"$%g$" % L)
        axes[1].plot(eta, G_rho, mk[L], fillstyle="none", ms=5, lw=1, c=clist[L])
    
    if rhoB == 0.03:
        axes[0].axvspan(0, 0.225, alpha=0.2, color="tab:blue", ec=None)
        axes[0].axvspan(0.225, 0.39, alpha=0.2, color="tab:orange", ec=None)
        axes[0].axvspan(0.39, 0.445, alpha=0.2, color="tab:blue", ec=None)
        axes[0].axvspan(0.445, 0.6, alpha=0.2, color="tab:grey", ec=None)
        axes[1].axvspan(0, 0.225, alpha=0.2, color="tab:blue", ec=None, label="bands")
        axes[1].axvspan(0.225, 0.39, alpha=0.2, color="tab:orange", ec=None, label="polar liquid")
        axes[1].axvspan(0.39, 0.445, alpha=0.2, color="tab:blue", ec=None)
        axes[1].axvspan(0.445, 0.6, alpha=0.2, color="tab:grey", ec=None, label="disordered gas")

if __name__ == "__main__":
    fig, axes = plt.subplots(2, 2, figsize=(8, 6), constrained_layout=True, sharey="row", sharex=True)

    dx = 40
    rhoB_arr = [0.03, 0.1]
    mk = {200:"-o", 400:"-s", 800: "->"}
    clist = {200: "tab:blue", 400: "tab:orange", 800: "tab:green"}
    for j, rhoB in enumerate(rhoB_arr):
        if rhoB == 0.03:
            L_arr = [200, 400, 800]
        else:
            L_arr = [200, 400, 800]
        for L in L_arr:
            eta, phi, G, G_rho, G_phi = read_eta_phi_G_Grho_Gphi(L, rhoB, dx=dx)
            if j == 0:
                label = None
            else:
                label=r"$L=%g$" % L
            axes[0, j].plot(eta, phi, mk[L], fillstyle="none", ms=5, label=label, lw=1, c=clist[L])
            axes[1, j].plot(eta, G_rho, mk[L], fillstyle="none", ms=5, lw=1, c=clist[L])
    
    dx, G_rho = get_order_para_varied_dx()
    ax_in = axes[1, 0].inset_axes([0.5, 0.1, 0.48, 0.4])
    ax_in.plot(dx, G_rho, "D", c="tab:red", ms=4, fillstyle="none")
    ax_in.axhline(2/3, linestyle="dashed", color="k")
    ax_in.text(0.75, 0.05, r"$\Delta x$", transform=ax_in.transAxes, fontsize="xx-large")
    ax_in.text(0.04, 0.78, r"$G_\rho$", transform=ax_in.transAxes, fontsize="xx-large")
    ax_in.set_title(r"$\eta=0.3, L=2000$", fontsize="x-large")
    ax_in.set_ylim(ymax=0.9)

    axes[0, 0].axvspan(0, 0.225, alpha=0.2, color="tab:blue", ec=None, label="bands")
    axes[0, 0].axvspan(0.225, 0.39, alpha=0.2, color="tab:orange", ec=None, label="polar liquid")
    axes[0, 0].axvspan(0.39, 0.445, alpha=0.2, color="tab:blue", ec=None)
    axes[0, 0].axvspan(0.445, 0.6, alpha=0.2, color="tab:grey", ec=None, label="disordered gas")

    axes[0, 0].legend(fontsize="xx-large", borderpad=0.2)

    axes[1, 0].axvspan(0, 0.225, alpha=0.2, color="tab:blue", ec=None)
    axes[1, 0].axvspan(0.225, 0.39, alpha=0.2, color="tab:orange", ec=None)
    axes[1, 0].axvspan(0.39, 0.445, alpha=0.2, color="tab:blue", ec=None)
    axes[1, 0].axvspan(0.445, 0.6, alpha=0.2, color="tab:grey", ec=None)

    axes[0, 1].axvspan(0, 0.425, alpha=0.2, color="tab:blue", ec=None)
    axes[0, 1].axvspan(0.425, 0.6, alpha=0.2, color="tab:grey", ec=None)
    axes[1, 1].axvspan(0, 0.425, alpha=0.2, color="tab:blue", ec=None)
    axes[1, 1].axvspan(0.425, 0.6, alpha=0.2, color="tab:grey", ec=None)

    axes[1, 0].set_ylim(ymin=-1, ymax=0.75)
    axes[0, 0].set_ylim(ymin=0, ymax=0.95)
    axes[0, 0].set_xlim(0.04, 0.51)
    axes[1, 0].axhline(2/3, linestyle="dashed", color="k")
    axes[1, 1].axhline(2/3, linestyle="dashed", color="k", label=r"$G_\rho=2/3$")
    axes[0, 1].legend(fontsize="xx-large")
    axes[0, 0].set_ylabel(r"$\langle \phi\rangle$", fontsize="xx-large")
    axes[1, 0].set_ylabel(r"$G_{\rho}$", fontsize="xx-large")
    axes[1, 0].set_xlabel(r"$\eta$", fontsize="xx-large")
    axes[1, 1].set_xlabel(r"$\eta$", fontsize="xx-large")
    axes[0, 0].set_title(r"$\bar{\rho}_B=0.03$", fontsize="xx-large")
    axes[0, 1].set_title(r"$\bar{\rho}_B=0.1$", fontsize="xx-large")

    axes[1, 1].legend(fontsize="xx-large", loc=(0.2, 0.7))

    axes[0, 0].text(0.02, 0.75, r"(a)", transform=axes[0,0].transAxes, fontsize="xx-large")
    axes[0, 1].text(0.02, 0.75, r"(b)", transform=axes[0,1].transAxes, fontsize="xx-large")
    axes[1, 0].text(0.02, 0.85, r"(c)", transform=axes[1,0].transAxes, fontsize="xx-large")
    axes[1, 0].text(0.02, 0.65, r"$\Delta x=40$", transform=axes[1,0].transAxes, fontsize="xx-large")
    # axes[1, 1].text(0.02, 0.65, r"$\Delta x=40$", transform=axes[1,1].transAxes, fontsize="xx-large")
    axes[1, 1].text(0.02, 0.85, r"(d)", transform=axes[1,1].transAxes, fontsize="xx-large")
    # plt.show()
    plt.savefig("fig/FIG2.pdf")
    plt.close()
