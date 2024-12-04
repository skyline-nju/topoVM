import numpy as np
import matplotlib.pyplot as plt
from add_line import add_line


root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"

def square_average(f, dx):
    if dx == 1:
        return np.mean(f), np.var(f)
    else:
        nrows, ncols = f.shape
        nrows_new = nrows // dx
        ncols_new = ncols // dx
        f_new = np.zeros((nrows_new, ncols_new))
        for j in range(nrows_new):
            for i in range(ncols_new):
                j0 = j * dx
                i0 = i * dx
                j1 = j0 + dx
                i1 = i0 + dx
                f_new[j, i] = np.sum(f[j0: j1, i0: i1])
        return np.mean(f_new), np.var(f_new)


def cal_GNF(L, rhoB=0.03, eta=0.3, dx=4):
    folder = f"{root_sohrab}/topoVM/dissenters/L{L}/coarse_grain_dx{dx}"
    basename = f"L{L}_{L}_d{rhoB:.04f}_e{eta:.3f}_r1_s3100.npz"
    with np.load(f"{folder}/{basename}", "r") as data:
        t, x, y, fields = data["t"], data["x"], data["y"], data["fields"]
    
    ncut = 10
    rho = fields[ncut:, 0]
    if dx == 4:
        if L == 2000:
            dn_arr = np.array([1, 2, 5, 10, 25, 50, 100, 125, 250])
        elif L == 1000:
            dn_arr = np.array([1, 2, 5, 10, 25, 50, 125])
    elif dx == 2:
        if L == 2000:
            dn_arr = np.array([1, 2, 4, 8, 10, 20, 25, 40, 50, 100, 125, 200, 250, 500])


    tot_frames = rho.shape[0]
    n_mean = np.zeros((tot_frames, dn_arr.size))
    n_var = np.zeros((tot_frames, dn_arr.size))

    for i_frame in range(tot_frames):
        for j, dn in enumerate(dn_arr):
            n_mean[i_frame, j], n_var[i_frame, j] = square_average(rho[i_frame] * dx**2, dn)
    
    # plt.loglog(np.mean(n_mean, axis=0), np.mean(n_var, axis=0), "o")
    # plt.show()
    # plt.close()

    fout = f"data/GNF/L{L}_d{rhoB:.2f}_e{eta:.2f}_r1_dx{dx}.npz"
    np.savez_compressed(fout, mean=np.mean(n_mean, axis=0), var=np.mean(n_var, axis=0))

def read_GNF(L, rhoB=0.03, eta=0.3, rho0=1, dx=4):
    fin = f"data/GNF/L{L}_d{rhoB:.2f}_e{eta:.2f}_r{rho0:g}_dx{dx}.npz"
    with np.load(fin, "r") as data:
        n_mean, n_var = data["mean"], data["var"]
    return n_mean, n_var


if __name__ == "__main__":
    # cal_GNF(2000, rhoB=0, eta=0.48)
    # cal_GNF(2000, rhoB=0, eta=0.5)
    # cal_GNF(2000, eta=0.45, dx=2, rhoB=0)
    # cal_GNF(1000, eta=0.3)
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    n_mean, n_var = read_GNF(2000, dx=2)
    ax.plot(n_mean, n_var, "o", fillstyle="none")
    n_mean, n_var = read_GNF(2000, dx=2, rhoB=0, eta=0.45)
    ax.plot(n_mean, n_var, "s", fillstyle="none")

    # n_mean, n_var = read_GNF(1000, eta=0.45)
    # ax.plot(n_mean, n_var, ">")

    # n_mean, n_var = read_GNF(2000, eta=0.50, rhoB=0)
    # ax.plot(n_mean, n_var, "v")

    ax.set_xscale("log")
    ax.set_yscale("log")

    add_line(ax, 0, 0, 1, 1.75)
    add_line(ax, 0, 0, 1, 1)

    plt.show()
    plt.close()

    




