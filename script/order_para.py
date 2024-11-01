import numpy as np
import matplotlib.pyplot as plt
import os

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"

def read_order_para_series(fin):
    with open(fin, "r") as f:
        lines = f.readlines()
        n = len(lines)-1
        # print("find", n, "lines in", os.path.basename(fin))
        t = np.zeros(n)
        theta = np.zeros(n)
        phi = np.zeros(n)
        for i, line in enumerate(lines[:-1]):
            s = line.rstrip("\n").split("\t")
            t[i], phi[i], theta[i] = float(s[0]), float(s[1]), float(s[2])
    return t, phi, theta

def cal_phi_G(phi_t, n_cut=2000):
    phi = phi_t[n_cut:]
    phi_m = np.mean(phi)
    phi2 = phi ** 2
    phi4 = phi2 ** 2
    G = 1 - np.mean(phi4) / (3 * np.mean(phi2) ** 2)
    return phi_m, G


def get_eta_arr(L, rhoB):
    prefix = f"{root_sohrab}/topoVM/dissenters"
    if L == 200:
        folder = f"{prefix}/L200_new"
        if rhoB == 0:
            eta_arr = np.array([
                0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5
            ])
        elif rhoB == 0.03:
            eta_arr = np.array([
                0, 0.05, 0.08, 0.1, 0.15, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
                0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43,
                0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5
            ])
        elif rhoB == 0.05:
            eta_arr = np.array([
                0, 0.05, 0.1, 0.15, 0.2, 0.25,
                0.3, 0.35, 0.40, 0.41, 0.42, 0.43,
                0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5
            ])
        elif rhoB == 0.1:
            eta_arr = np.array([
                0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5
            ])
    elif L == 400:
        folder = f"{prefix}/L400_new"
        if rhoB == 0:
            eta_arr = np.array([
                0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
                0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5
            ])
        elif rhoB == 0.03:
            eta_arr = np.array([
                0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
                0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
                0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
                0.3, 0.35, 0.38,
                0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5
            ])
        elif rhoB == 0.05:
            eta_arr = np.array([
                0, 0.01, 0.02, 0.03, 0.04, 0.05, 
                0.1, 0.15,
                0.2, 0.25,
                0.3, 0.35, 0.38,
                0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5
            ])
        elif rhoB == 0.1:
            eta_arr = np.array([
                0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5
            ])
    elif L == 800:
        folder = f"{prefix}/L800_new2"
        if rhoB == 0.03:
            eta_arr = np.array([
                0.05, 0.06, 0.07, 0.08, 0.09, 
                0.1, 0.15,
                0.2, 0.25,
                0.3, 0.35, 0.38,
                0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5
            ])
        elif rhoB == 0.1:
            eta_arr = np.array([
                0.05,
                0.1, 0.15,
                0.2, 0.25,
                0.3, 0.35,
                0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49
            ])
        elif rhoB == 0:
            eta_arr = np.array([
                0.05,
                0.1, 0.15,
                0.2, 0.25,
                0.3, 0.35,
                0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49
            ])
    elif L == 1600:
        folder = f"{prefix}/L1600"
        if rhoB == 0.03:
            eta_arr = np.array([0.2, 0.25, 0.3, 0.35, 0.4, 0.42])
    return eta_arr, folder
        

def get_phi_G_arr(L, rhoB, rho0=1., seed=3100, ncut=2000):
    eta_arr, folder = get_eta_arr(L, rhoB)
    if L == 200 or L == 400 or  L == 1600:
        folder = f"{folder}/phi"
    else:
        folder = f"{folder}/op"
    n = eta_arr.size
    phi_arr, G_arr = np.zeros((2, n))

    for i, eta in enumerate(eta_arr):
        if L == 200 or L == 400 or L == 1600:
            fin = f"{folder}/{L:d}_{L:d}_{rhoB:.4f}_{eta:.3f}_{rho0:.3f}_{seed:d}.dat"
        elif L == 800:
            t0 = 0
            fin = f"{folder}/L{L:d}_{L:d}_d{rhoB:.4f}_e{eta:.3f}_r{rho0:g}_s{seed:d}_t{t0:08d}.dat"
        t, phi, theta = read_order_para_series(fin)
        phi_arr[i], G_arr[i] = cal_phi_G(phi, ncut)
    return eta_arr, phi_arr, G_arr


def get_Grho_Gphi_arr(L, bar_rhoB, dx=25, rho0=1., seed=3100):
    def cal_G(f):
        f4_m = np.mean(f ** 4)
        f2_m = np.mean(f ** 2)
        G = 1 - f4_m / (3 * f2_m ** 2)
        return G

    eta_arr, folder = get_eta_arr(L, bar_rhoB)
    G_rho_arr = np.zeros_like(eta_arr)
    G_phi_arr = np.zeros_like(eta_arr)
    folder = f"{folder}/coarse_grain_dx{dx}"
    for i, eta in enumerate(eta_arr):
        basename = f"L{L:d}_{L:d}_d{bar_rhoB:.4f}_e{eta:.3f}_r{rho0:g}_s{seed:d}.npz"
        fname = f"{folder}/{basename}"
        with np.load(fname, "r") as data:
            t = data["t"]
            x = data["x"]
            fields = data["fields"]

            if t.size < 75:
                ncut = 10
            else:
                ncut = 20
            rhoA = fields[ncut:, 0]
            rhoB = fields[ncut:, 1]
            mA_x = fields[ncut:, 2]
            mA_y = fields[ncut:, 4]

        G_rho_arr[i] = cal_G(rhoA)
        phi = np.sqrt(mA_x**2 + mA_y**2) / rhoA
        G_phi_arr[i] = cal_G(phi)


    return eta_arr, G_rho_arr, G_phi_arr


def cal_order_para_all(dx):
    rhoB_arr = [0, 0.03, 0.05, 0.1]
    for rhoB in rhoB_arr:
        if rhoB == 0.03:
            if dx == 40:
                L_arr = [200, 400, 800, 1600]
            else:
                L_arr = [200, 400, 800]
        elif rhoB == 0.1 or rhoB == 0.:
            L_arr = [200, 400, 800]
        else:
            L_arr = [200, 400]
        for L in L_arr:
            eta, phi, G = get_phi_G_arr(L, rhoB)
            eta, G_rho, G_phi = get_Grho_Gphi_arr(L, rhoB, dx=dx)
            fout = f"data/order_para/L{L:g}_rhoB{rhoB:g}_dx{dx:g}.npz"
            np.savez_compressed(fout, eta=eta, phi=phi, G=G, G_rho=G_rho, G_phi=G_phi)
    

def read_eta_phi_G_Grho_Gphi(L, rhoB, dx):
    fin = f"data/order_para/L{L:g}_rhoB{rhoB:g}_dx{dx:g}.npz"
    with np.load(fin, "r") as data:
        eta = data["eta"]
        phi = data["phi"]
        G = data["G"]
        G_rho = data["G_rho"]
        G_phi = data["G_phi"]
    return eta, phi, G, G_rho, G_phi


def plot_phi_G_G_rho():
    fig, axes = plt.subplots(3, 3, sharex=True, constrained_layout=True)
    rhoB_arr = [0, 0.03, 0.05]

    for j, rhoB in enumerate(rhoB_arr):
        if rhoB == 0.03:
            L_arr = [200, 400, 800]
        else:
            L_arr = [200, 400]
        for L in L_arr:
            eta, phi, G, G_rho, G_phi = read_eta_phi_G_Grho_Gphi(L, rhoB)
            axes[j, 0].plot(eta, phi, "-o")
            axes[j, 1].plot(eta, G, '-o')
            axes[j, 2].plot(eta, G_rho, '-o')
    plt.show()
    plt.close()


if __name__ == "__main__":
    dx = 40
    cal_order_para_all(dx=dx)

    # fig, axes = plt.subplots(2, 3, sharex=True, constrained_layout=True, sharey="row")

    # rhoB_arr = [0., 0.03, 0.1]
    # for j, rhoB in enumerate(rhoB_arr):
    #     if rhoB == 0.03:
    #         L_arr = [200, 400, 800]
    #     else:
    #         L_arr = [200, 400]
    #     for L in L_arr:
    #         eta, phi, G, G_rho, G_phi = read_eta_phi_G_Grho_Gphi(L, rhoB, dx=dx)
    #         axes[0, j].plot(eta, phi, "-o")
    #         axes[1, j].plot(eta, G_rho, '-o')
    #         # axes[1, j].plot(eta, G_phi, '-s')
    # axes[1, 0].set_ylim(-1, 1)
    # plt.show()
    # plt.close()
