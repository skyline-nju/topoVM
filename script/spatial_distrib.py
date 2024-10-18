import numpy as np
import matplotlib.pyplot as plt

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def show_distrib(fname, ncut=20, bins=100):
    with np.load(fname, "r") as data:
        t = data["t"]
        x = data["x"]
        fields = data["fields"]

        rhoA = fields[ncut:, 0]
        rhoB = fields[ncut:, 1]

        print(fields.shape)

        hist_A, bin_edges = np.histogram(rhoA, bins=bins, density=True)
        rA = (bin_edges[1:] + bin_edges[:-1]) * 0.5

        hist_B, bin_edges = np.histogram(rhoB, bins=bins, density=True)
        rB = (bin_edges[1:] + bin_edges[:-1]) * 0.5

        plt.plot(rA, hist_A, "-o")
        plt.plot(rB, hist_B, "-s")

        plt.show()
        plt.close()


def cal_G_rho(pB, dx=25, ncut=20, L=400):
    if pB == 0.05:
        eta_arr = np.array([0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5])
    elif pB == 0.03:
        eta_arr = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45,
                            0.5])
    G_rho_arr = np.zeros_like(eta_arr)

    folder = f"{root_sohrab}/topoVM/dissenters/L{L:d}_new/coarse_grain_dx{dx}"
    for i, eta in enumerate(eta_arr):
        basename = f"L400_400_d{pB:.4f}_e{eta:.3f}_r1_s3100.npz"
        fname = f"{folder}/{basename}"
        with np.load(fname, "r") as data:
            t = data["t"]
            x = data["x"]
            fields = data["fields"]

            rhoA = fields[ncut:, 0]
            rhoB = fields[ncut:, 1]
        
        rhoA4_m = np.mean(rhoA ** 4)
        rhoA2_m = np.mean(rhoA ** 2)
        G_rho_arr[i] = 1 - rhoA4_m / (3 * rhoA2_m ** 2)
    return eta_arr, G_rho_arr


def cal_G_phi(pB, dx=25, ncut=20, L=400):
    if pB == 0.05:
        eta_arr = np.array([0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5])
    elif pB == 0.03:
        eta_arr = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45,
                            0.5])
    G_phi_arr = np.zeros_like(eta_arr)

    folder = f"{root_sohrab}/topoVM/dissenters/L{L:d}_new/coarse_grain_dx{dx}"
    for i, eta in enumerate(eta_arr):
        basename = f"L400_400_d{pB:.4f}_e{eta:.3f}_r1_s3100.npz"
        fname = f"{folder}/{basename}"
        with np.load(fname, "r") as data:
            t = data["t"]
            x = data["x"]
            fields = data["fields"]

            rhoA = fields[ncut:, 0]
            rhoB = fields[ncut:, 1]
            mA_x = fields[ncut:, 2]
            mA_y = fields[ncut:, 4]

        # phi = np.sqrt(mA_x ** 2 + mA_y ** 2)
        mask = rhoA > 0
        mA_x[mask] /= rhoA[mask]
        mA_y[mask] /= rhoA[mask]

        phi = np.arctan2(mA_y, mA_x)
        
        phiA4_m = np.mean(phi ** 4)
        phiA2_m = np.mean(phi ** 2)
        G_phi_arr[i] = 1 - phiA4_m / (3 * phiA2_m ** 2)
    return eta_arr, G_phi_arr


if __name__ == "__main__":
    # dx = 25
    # pB = 0.05

    # ncut = 25

    # eta_arr, G_rho = cal_G_rho(0.03)
    # plt.plot(eta_arr, G_rho, "-o")
    # eta_arr, G_rho = cal_G_rho(0.05)
    # plt.plot(eta_arr, G_rho, "-s")

    # plt.show()
    # plt.close()

    dx=25

    folder = f"{root_sohrab}/topoVM/dissenters/L2000/coarse_grain_dx{dx:d}"

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    eta_arr = np.array([0.1, 0.3, 0.42])
    ncut_arr = [6, 40, 30]
    for i, eta in enumerate(eta_arr):
        fname = f"{folder}/L2000_2000_d0.0300_e{eta:.3f}_r1_s3100.npz"
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
        ax1.plot(rho_arr, pdf_rho, "-o")

        phi_arr = np.zeros_like(rho_arr)
        for i in range(rho_arr.size):
            mask = np.logical_and(rhoA >= bin_edges[i], rhoA < bin_edges[i+1])
            if np.sum(mask) > 0:
                phi_arr[i] = np.mean(phi[mask])
        
        ax3.plot(rho_arr, phi_arr, "-o")

        bins = 100
        pdf_phi, bin_edges = np.histogram(phi, bins=bins, density=True)
        phi_arr = (bin_edges[1:] + bin_edges[:-1]) * 0.5
        ax2.plot(phi_arr, pdf_phi, "-o")

        # ax3.plot()
   
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    # ax3.set_xscale("log")
    # ax3.set_yscale("log")
    plt.show()
    plt.close()



