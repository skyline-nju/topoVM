import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy import linalg

def get_J_kl(k: int, l: int)-> float:
    if 2 * l - k == 0:
        return 1.
    else:
        if k % 2 == 0:
            return 0.
        elif (2 * l - k ) % 4 == 1:
            return 1 / (np.pi * (l - k/2))
        else:
            return -1 / (np.pi * (l - k/2))


def get_P_k(k: int, eta: float)-> float:
    return np.exp(-0.5 * (k * eta)**2)


def get_eta_c(rhoA, rhoB, alpha=1.):
    rho0 = rhoA + rhoB
    phiB = rhoB / rho0
    if not isinstance(phiB, np.ndarray):
        phiB = np.array([phiB])
        is_scalar = True
    else:
        is_scalar = False
    mask = 2  * phiB - (4 - np.pi) < 0
    eta = np.zeros_like(phiB)
    eta[mask] = np.sqrt(
        2 * np.log ((np.pi + 4 * alpha - 2 * alpha * phiB[mask])/(np.pi * (1+alpha)))
    )
    if is_scalar:
        eta = eta[0]
    return eta


def func_HS(x, rhoA, rhoB, eta, alpha):
    K = x.size - 1
    y = np.zeros(K+1)
    rho0 = rhoA + rhoB
    for k in range(0, y.size):
        P_k = get_P_k(k, eta)
        y[k] = (P_k - 1 - alpha) * x[k]
        S = 0.
        for l in range(-K, K+1):
            m = abs(k - l)
            if m <= K:
                S += get_J_kl(k, l) * x[abs(l)] * x[m]
        y[k] += alpha / rho0 * P_k * (S + get_J_kl(k, 0) * x[k] * rhoB)
    return y

def get_f_bar(K, eta, rhoA, rhoB, alpha):
    x0 = np.zeros(K+1)
    x0[0] = rhoA
    for i in range(1, K+1):
        x0[i] = 0.8 ** i
    f_bar = fsolve(func_HS, x0, args=(rhoA, rhoB, eta, alpha))
    return f_bar

def get_M_cross_base(K: int):
    M_gg_base = np.zeros((K+1, K+1), complex)
    M_gh_base = np.zeros((K+1, K+1), complex)
    M_hg_base = np.zeros((K+1, K+1), complex)
    M_hh_base = np.zeros((K+1, K+1), complex)

    # k=0 row
    M_gg_base[0, 1] = 2
    M_gh_base[0, 1] = 2
    # M_hg_q[0, 1] = 0.
    # M_hh_q[0, 1] = 0.

    # k=1, 2, ..., K-1
    for k in range(1, K):
        M_gg_base[k, k+1] = 1
        M_gg_base[k, k-1] = 1
        M_gh_base[k, k+1] = 1
        M_gh_base[k, k-1] = -1

        M_hg_base[k, k+1] = -1
        M_hg_base[k, k-1] = 1
        M_hh_base[k, k+1] = 1
        M_hh_base[k, k-1] = 1

    # k=K
    M_gg_base[K, K-1] = 1
    M_gh_base[K, K-1] = -1
    M_hg_base[K, K-1] = 1
    M_hh_base[K, K-1] = 1
    return M_gg_base, M_gh_base, M_hg_base, M_hh_base


def get_M_q(qx, qy, M_gg_base, M_gh_base, M_hg_base, M_hh_base):
    M_gg_q = -0.5j * qx * M_gg_base
    M_gh_q = -0.5j * qy * M_gh_base
    M_hg_q = -0.5j * qy * M_hg_base
    M_hh_q = -0.5j * qx * M_hh_base
    return M_gg_q, M_gh_q, M_hg_q, M_hh_q


def get_M_0(eta, rhoA, rhoB, alpha, K):
    rho0 = rhoA + rhoB
    f_bar = get_f_bar(K, eta, rhoA, rhoB, alpha)
    Pk = np.array([get_P_k(k, eta) for k in range(K+1)])

    diag_arr = Pk - 1 - alpha
    M_gg_0 = np.diag(diag_arr)
    M_gg_0[:, 0] += diag_arr * f_bar / rho0

    M_hh_0 = np.copy(M_gg_0)

    for k in range(0, K+1):
        tmp = alpha / rho0 * Pk[k]
        # j = 0 column
        M_k0 = tmp * (get_J_kl(k, 0) + get_J_kl(k, k)) * f_bar[k]
        M_gg_0[k, 0] += M_k0
        M_hh_0[k, 0] += M_k0
        for j in range(1, K+1):
            m1 = abs(k-j)
            sum1 = (get_J_kl(k, j) + get_J_kl(k, k-j)) * f_bar[m1]

            l = -j
            m2 = abs(k-l)
            if m2 <= K:
                sum2 = (get_J_kl(k, l) + get_J_kl(k, k-l)) * f_bar[m2]
            else:
                sum2 = 0
            
            M_gg_0[k, j] += tmp * (sum1 + sum2)
            M_hh_0[k, j] += tmp * (sum1 - sum2)
        M_gg_0[k, k] += tmp * get_J_kl(k, 0) * rhoB
        M_hh_0[k, k] += tmp * get_J_kl(k, 0) * rhoB
    return M_gg_0, M_hh_0, f_bar


def assemble_matrix(M_gg_0, M_hh_0, M_gg_base, M_gh_base, M_hg_base, M_hh_base, qx, qy):
    M_gg_q, M_gh_q, M_hg_q, M_hh_q = get_M_q(qx, qy, M_gg_base, M_gh_base, M_hg_base, M_hh_base)
    n = M_gg_0.shape[0]
    M = np.zeros((2 * n, 2 * n), complex)
    M[:n, :n] = M_gg_0 + M_gg_q
    M[:n, n:] = M_gh_q
    M[n:, :n] = M_hg_q
    M[n:, n:] = M_hh_0 + M_hh_q
    return M

def get_M(eta, rhoA, rhoB, alpha, K, qx, qy):
    M_gg_0, M_hh_0, f_bar = get_M_0(eta, rhoA, rhoB, alpha, K)
    M_gg_base, M_gh_base, M_hg_base, M_hh_base = get_M_cross_base(K)
    M = assemble_matrix(M_gg_0, M_hh_0, M_gg_base, M_gh_base, M_hg_base, M_hh_base, qx, qy)
    return M, f_bar


def plot_sigma_vs_rhoB():
    rho0 = 1
    rhoB = 0.08
    rhoA = rho0-rhoB
    alpha = 1
    qx = 1e-3
    qy = 0

    fig, (ax1, ax2) = plt.subplots(2, 1,sharex=True, figsize=(4, 8), constrained_layout=True)
    ax1.axhline(0, linestyle="dashed", c="tab:red")

    K = 20
    eta_arr = np.linspace(0, 0.55, 100)
    sigma = np.zeros_like(eta_arr)
    f1_arr = np.zeros_like(eta_arr)

    for rhoB in [0.05, 0.1, 0.12, 0.13, 0.14, 0.15]:
        for i, eta in enumerate(eta_arr):
            M, f_bar = get_M(eta, rhoA, rhoB, alpha, K, qx, qy)
            sigma_re = linalg.eigvals(M).real
            sigma[i] = np.max(sigma_re)
            f1_arr[i] = f_bar[1]
        ax1.plot(eta_arr, sigma, "-")
        ax2.plot(eta_arr, f1_arr, "-")
        print(sigma.max(),sigma.min())
    ax1.set_ylim(0, 1e-7)
    plt.show()
    plt.close()


def cal_sigma_qx_fixed_rho0(K=50, n=50):
    rhoB_arr = np.linspace(0, 0.5, n)
    eta_arr = np.linspace(0, 0.5, n)

    qx = 1e-3
    qy = 0

    rho0 = 1
    alpha = 1
    
    sigma_arr = np.zeros((eta_arr.size, rhoB_arr.size))

    for j, eta in enumerate(eta_arr):
        print("j=", j)
        for i, rhoB in enumerate(rhoB_arr):
            rhoA = rho0 - rhoB
            M, f_bar = get_M(eta, rhoA, rhoB, alpha, K, qx, qy)
            sigma_re = linalg.eigvals(M).real
            sigma_arr[j, i] = np.max(sigma_re)
    
    f_out = f"data/r{rho0:g}_qx{qx:.4f}_K{K:d}_n{n:d}.npz"
    np.savez_compressed(f_out, rhoB=rhoB_arr, eta=eta_arr, sigma=sigma_arr)


def cal_sigma_qx_fixed_rhoB(K=50, n=50, skip_neg_mu=False):
    rhoA_arr = np.linspace(0, 1.2, n)
    eta_arr = np.linspace(0, 0.5, n)

    qx = 1e-3
    qy = 0

    rhoB = 0.03

    alpha = 1
    
    sigma_arr = np.zeros((eta_arr.size, rhoA_arr.size))

    for j, eta in enumerate(eta_arr):
        print("j=", j)
        for i, rhoA in enumerate(rhoA_arr):
            rho0 = rhoA + rhoB
            eta_c = get_eta_c(rhoA, rhoB, alpha=alpha)
            if eta > eta_c and skip_neg_mu:
                sigma_arr[j, i] = 0
            else:
                M, f_bar = get_M(eta, rhoA, rhoB, alpha, K, qx, qy)
                sigma_re = linalg.eigvals(M).real
                sigma_arr[j, i] = np.max(sigma_re)
    
    f_out = f"data/rB{rhoB:g}_qx{qx:.4f}_K{K:d}_n{n:d}.npz"
    np.savez_compressed(f_out, rhoA=rhoA_arr, eta=eta_arr, sigma=sigma_arr)


def PD_rhoB_eta():
    # cal_sigma_qx_fixed_rho0(K=20, n=200)
    # cal_sigma_qx_fixed_rho0(K=50, n=200)
    # cal_sigma_qx_fixed_rho0(K=100, n=200)


    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True, sharex=True, sharey=True)
    K_arr = [20, 50, 100]
    n = [200, 200, 200]
    for i, K in enumerate(K_arr):
        if n[i] == 50:
            fin = f"data/r1_qx0.0010_K{K:d}.npz"
        else:
            fin = f"data/r1_qx0.0010_K{K:d}_n{n[i]:d}.npz"
        with np.load(fin, "r") as data:
            rhoB = data["rhoB"]
            eta = data["eta"]
            sigma = data["sigma"]
            state = np.zeros_like(sigma)
            mask = sigma > 1e-9
            state[mask] = 1
    
        dx = rhoB[1] - rhoB[0]
        dy = eta[1] - eta[0]
        extent = [rhoB[0]-dx/2, rhoB[-1]+dx/2, eta[0]-dy/2, eta[-1]+dy/2]
        axes[i].imshow(state, origin="lower",extent=extent)


        eta_c = get_eta_c(1-rhoB, rhoB, alpha=1)
        axes[i].plot(rhoB, eta_c, c="r")
        axes[i].set_xlim(rhoB[0], rhoB[-1])
        axes[i].set_ylim(eta[0], eta[-1])

    plt.show()
    plt.close()


def plot_rhoB_eta_panel(ax, K, n, rho0=1):
    if n == 50:
        fin = f"data/r1_qx0.0010_K{K:d}.npz"
    else:
        fin = f"data/r1_qx0.0010_K{K:d}_n{n:d}.npz"
    with np.load(fin, "r") as data:
        rhoB = data["rhoB"]
        eta = data["eta"]
        sigma = data["sigma"]
        state = np.zeros_like(sigma)
        mask = sigma > 1e-9
        state[mask] = 1

    dx = rhoB[1] - rhoB[0]
    dy = eta[1] - eta[0]
    extent = [rhoB[0]-dx/2, rhoB[-1]+dx/2, eta[0]-dy/2, eta[-1]+dy/2]
    ax.imshow(state, origin="lower",extent=extent, aspect="auto")
    eta_c = get_eta_c(1-rhoB, rhoB, alpha=1)
    ax.plot(rhoB, eta_c, c="r")
    ax.set_xlim(rhoB[0], rhoB[-1])
    ax.set_ylim(eta[0], eta[-1])

def plot_rhoA_eta_panel(ax, K, n, rhoB=0.03):
    fin = f"data/rB{rhoB:g}_qx0.0010_K{K:d}_n{n:d}.npz"
    with np.load(fin, "r") as data:
        eta = data["eta"]
        rhoA = data["rhoA"]
        sigma = data["sigma"]
        state = np.zeros_like(sigma)
        mask = sigma > 1e-9
        state[mask] = 1
        dx = rhoA[1] - rhoA[0]
        dy = eta[1] - eta[0]
        extent = [rhoA[0]-dx/2, rhoA[-1]+dx/2, eta[0]-dy/2, eta[-1]+dy/2]
        ax.imshow(state, origin="lower", extent=extent, aspect="auto")
        eta_c = get_eta_c(rhoA, rhoB, alpha=1)
        ax.plot(rhoA, eta_c, c="tab:red")

def PD_rhoA_eta():
    # cal_sigma_qx_fixed_rhoB(K=20, n=200, skip_neg_mu=True)
    # cal_sigma_qx_fixed_rhoB(K=50, n=200, skip_neg_mu=True)
    # cal_sigma_qx_fixed_rhoB(K=100, n=200, skip_neg_mu=True)

    fig, axes = plt.subplots(1, 3, constrained_layout=True, figsize=(12, 4), sharex=True, sharey=True)

    rhoB = 0.03
    K_arr = [20, 50, 100]
    n_arr = [200, 200, 200]
    for i, ax in enumerate(axes):
        plot_rhoA_eta_panel(ax, K_arr[i], n_arr[i])


    plt.show()
    plt.close()


if __name__ == "__main__":
    PD_rhoA_eta()
    # PD_rhoB_eta()