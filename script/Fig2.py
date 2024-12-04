import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import image as mpimg
from order_para import read_order_para_series
from GNF import read_GNF
from band_profile import plot_band_profiles_varied_rhoB, plot_band_profiles_varied_eta
from add_line import add_line

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


def get_phi_arr():
    L_arr = np.array([200, 400, 800, 1000, 1600, 2000])
    ncut = 2000
    phi_arr = np.zeros(L_arr.size)
    for i, L in enumerate(L_arr):
        fin = f"data/order_para/{L}_{L}_0.0300_0.300_1.000_3100.dat"
        t, phi, theta = read_order_para_series(fin)
        phi_arr[i] = np.mean(phi[ncut:])
    return L_arr, phi_arr


def one_column():
    fig, (ax2, ax3, ax1) = plt.subplots(1, 3, figsize=(8, 3.4), layout="constrained")


    n_mean, n_var = read_GNF(2000, dx=2)
    ax1.plot(n_mean, n_var, "-o", fillstyle="none", c="tab:orange", ms=4)
    n_mean, n_var = read_GNF(2000, dx=2, rhoB=0, eta=0.44)
    ax1.plot(n_mean, n_var, "-D", fillstyle="none", c="tab:green", ms=4)
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    dx=40
    folder = f"{root_sohrab}/topoVM/dissenters/L2000/coarse_grain_dx{dx:d}"

    eta_arr = [0.1, 0.3, 0.42]
    ncut_arr = [6, 40, 30, 20]

    fnames = [f"{folder}/L2000_2000_d0.0300_e{eta:.3f}_r1_s3100.npz" for eta in eta_arr]
    fnames.append(f"{folder}/L2000_2000_d0.0000_e0.440_r1_s3100.npz")
    eta_arr.append(0.44)
    clist = ["tab:blue", "tab:orange", "tab:cyan", "tab:green"]
    mk_list = ["-s", "-o", "->", "-D"]
    label_list = [r"$(\bar{\rho}_B,\eta)=(0.03, 0.1)$", r"$(0.03, 0.3)$", r"$(0.03, 0.42)$", r"$(0, 0.44)$"]
    for i, eta in enumerate(eta_arr):
        fname = fnames[i]
        ncut = ncut_arr[i]
        with np.load(fname, "r") as data:
            t = data["t"]
            x = data["x"]
            fields = data["fields"]

            rhoA = fields[ncut:, 0]
            rhoB = fields[ncut:, 1]
            mA_x = fields[ncut:, 2]
            mA_y = fields[ncut:, 4]
            phi = np.zeros_like(rhoA)
            mask = rhoA > 0.
            phi[mask] = np.sqrt(mA_x[mask]**2 + mA_y[mask]**2) / rhoA[mask]

        bins = 100
        pdf_rho, bin_edges = np.histogram(rhoA, bins=bins, density=True)
        rho_arr = (bin_edges[1:] + bin_edges[:-1]) * 0.5
        ax2.plot(rho_arr, pdf_rho, mk_list[i], fillstyle="none", c=clist[i], ms=4)

        phi_arr = np.zeros_like(rho_arr)
        for j in range(rho_arr.size):
            mask = np.logical_and(rhoA >= bin_edges[j], rhoA < bin_edges[j+1])
            if np.sum(mask) > 0:
                phi_arr[j] = np.mean(phi[mask])
        if i == 3:
            rho_arr = rho_arr[2:]
            phi_arr = phi_arr[2:]
        ax3.plot(rho_arr, phi_arr, mk_list[i], fillstyle="none", c=clist[i], ms=4, label=label_list[i])

    ax3.axvline(1, linestyle="dashed", c="tab:grey")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax2.set_xlim(xmin=6e-2, xmax=10)
    ax2.set_ylim(1e-3)
    ax3.set_xlim(0, 1.55)
    ax3.set_ylim(0.3, 0.95)

    # ax1.legend(title=r"$(\bar{\rho}_B,\eta)=$", fontsize="large", frameon=False, title_fontsize="large")
    # ax3.legend(title=r"$(\bar{\rho}_B,\eta)=$", fontsize="large", frameon=False, title_fontsize="large")

    ax1.set_xlabel(r"$\langle n\rangle$", fontsize="xx-large")
    ax1.set_ylabel(r"$\Delta n$", fontsize="xx-large")
    ax2.set_ylabel(r"$P(\rho_A)$", fontsize="xx-large")
    ax2.set_xlabel(r"$\rho_A(\mathbf{r})$", fontsize="xx-large")
    ax3.set_ylabel(r"$\phi(\mathbf{r})$", fontsize="xx-large")
    ax3.set_xlabel(r"$\rho_A(\mathbf{r})$", fontsize="xx-large")

    ax1.set_title(r"(c)", fontsize="xx-large")
    ax2.set_title(r"(a)", fontsize="xx-large")
    ax3.set_title(r"(b)", fontsize="xx-large")

    fig.legend(loc='outside lower center', ncols=4, fontsize="x-large", columnspacing=1.5, handletextpad=0.4)
    add_line(ax1, 0, 0.05, 0.9, 1.75, "1.75", xl=0.4, yl=0.65)
    add_line(ax2, 0, 0.95, 0.99, -1, "-1", xl=0.42, yl=0.5)
    plt.show()
    # plt.savefig("fig/FIG3.pdf")
    plt.close()

if __name__ == "__main__":
    fig, (ax2, ax3, ax1) = plt.subplots(1, 3, figsize=(8, 3.4), layout="constrained")
    n_mean, n_var = read_GNF(2000, dx=2)
    ax1.plot(n_mean, n_var, "-o", fillstyle="none", c="tab:orange", ms=4)
    n_mean, n_var = read_GNF(2000, dx=2, rhoB=0, eta=0.45)
    ax1.plot(n_mean, n_var, "-D", fillstyle="none", c="tab:green", ms=4)
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    dx=40
    folder = f"{root_sohrab}/topoVM/dissenters/L2000/coarse_grain_dx{dx:d}"

    eta_arr = [0.1, 0.3, 0.42]
    ncut_arr = [6, 40, 30, 20]

    fnames = [f"{folder}/L2000_2000_d0.0300_e{eta:.3f}_r1_s3100.npz" for eta in eta_arr]
    fnames.append(f"{folder}/L2000_2000_d0.0000_e0.450_r1_s3100.npz")
    eta_arr.append(0.45)
    clist = ["tab:blue", "tab:orange", "tab:cyan", "tab:green"]
    mk_list = ["-s", "-o", "->", "-D"]
    label_list = [r"$(\bar{\rho}_B,\eta)=(0.03, 0.1)$", r"$(0.03, 0.3)$", r"$(0.03, 0.42)$", r"$(0, 0.45)$"]
    for i, eta in enumerate(eta_arr):
        fname = fnames[i]
        ncut = ncut_arr[i]
        with np.load(fname, "r") as data:
            t = data["t"]
            x = data["x"]
            fields = data["fields"]

            rhoA = fields[ncut:, 0]
            rhoB = fields[ncut:, 1]
            mA_x = fields[ncut:, 2]
            mA_y = fields[ncut:, 4]
            phi = np.zeros_like(rhoA)
            mask = rhoA > 0.
            phi[mask] = np.sqrt(mA_x[mask]**2 + mA_y[mask]**2) / rhoA[mask]

        bins = 100
        pdf_rho, bin_edges = np.histogram(rhoA, bins=bins, density=True)
        rho_arr = (bin_edges[1:] + bin_edges[:-1]) * 0.5
        ax2.plot(rho_arr, pdf_rho, mk_list[i], fillstyle="none", c=clist[i], ms=4)

        phi_arr = np.zeros_like(rho_arr)
        for j in range(rho_arr.size):
            mask = np.logical_and(rhoA >= bin_edges[j], rhoA < bin_edges[j+1])
            if np.sum(mask) > 0:
                phi_arr[j] = np.mean(phi[mask])
        if i == 3:
            rho_arr = rho_arr[2:]
            phi_arr = phi_arr[2:]
        ax3.plot(rho_arr, phi_arr, mk_list[i], fillstyle="none", c=clist[i], ms=4, label=label_list[i])

    ax3.axvline(1, linestyle="dashed", c="tab:grey")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax2.set_xlim(xmin=5e-2, xmax=10)
    ax2.set_ylim(1e-3, 2)
    ax3.set_xlim(0, 1.55)
    ax3.set_ylim(0.3, 0.95)

    # ax1.legend(title=r"$(\bar{\rho}_B,\eta)=$", fontsize="large", frameon=False, title_fontsize="large")
    # ax3.legend(title=r"$(\bar{\rho}_B,\eta)=$", fontsize="large", frameon=False, title_fontsize="large")

    # ax1.set_xlabel(r"$\langle n\rangle$", fontsize="xx-large")
    # ax1.set_ylabel(r"$\Delta n$", fontsize="xx-large")
    # ax2.set_ylabel(r"$P(\rho_A)$", fontsize="xx-large")
    # ax2.set_xlabel(r"$\rho_A(\mathbf{r})$", fontsize="xx-large")
    # ax3.set_ylabel(r"$\phi(\mathbf{r})$", fontsize="xx-large")
    # ax3.set_xlabel(r"$\rho_A(\mathbf{r})$", fontsize="xx-large")

    ax2.text(0.02, 0.75, r"$P(\rho_A)$", fontsize="xx-large", transform=ax2.transAxes, rotation=90)
    ax3.text(0.03, 0.75, r"$\phi(\mathbf{r})$", fontsize="xx-large", transform=ax3.transAxes, rotation=90)
    ax1.text(0.03, 0.8, r"$\Delta n$", fontsize="xx-large", transform=ax1.transAxes, rotation=90)

    ax2.text(0.82, 0.04, r"$\rho_A$", fontsize="xx-large", transform=ax2.transAxes)
    ax3.text(0.7, 0.04, r"$\rho_A(\mathbf{r})$", fontsize="xx-large", transform=ax3.transAxes)
    ax1.text(0.77, 0.04, r"$\langle n\rangle$", fontsize="xx-large", transform=ax1.transAxes)

    ax1.set_title(r"(c)", fontsize="xx-large")
    ax2.set_title(r"(a)", fontsize="xx-large")
    ax3.set_title(r"(b)", fontsize="xx-large")

    fig.legend(loc='outside lower center', ncols=4, fontsize="x-large", columnspacing=1.5, handletextpad=0.4)
    add_line(ax1, 0., 0.05, 0.9, 1.75, "1.75", xl=0.4, yl=0.65)
    add_line(ax2, 0.15, 0.9, 0.99, -1, "-1", xl=0.42, yl=0.5)
    plt.show()
    # plt.savefig("fig/FIG2.pdf")
    plt.close()
