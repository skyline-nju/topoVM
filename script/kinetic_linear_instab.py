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

def get_f_bar(K, eta, rhoA, rhoB, alpha, f_bar_pre=None):
    if f_bar_pre is None:
        x0 = np.zeros(K+1)
        x0[0] = rhoA
        for i in range(1, K+1):
            x0[i] = 0.8 ** i
    else:
        x0 = f_bar_pre
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


def get_M_0(eta, rhoA, rhoB, alpha, K, f_bar_pre=None):
    rho0 = rhoA + rhoB
    f_bar = get_f_bar(K, eta, rhoA, rhoB, alpha, f_bar_pre=f_bar_pre)

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


def cal_sigma_q_theta_fixed_rho0(K=50, n=50, n_theta=30, skip_neg_mu=False):
    rhoB_arr = np.linspace(0, 0.5, n)
    eta_arr = np.linspace(0, 0.5, n)
    q = 1e-3
    rho0 = 1
    alpha = 1
    sigma_arr = np.zeros((eta_arr.size, rhoB_arr.size))
    theta_arr = np.zeros_like(sigma_arr)

    theta_1D = np.linspace(0, np.pi/2, n_theta)

    for j, eta in enumerate(eta_arr):
        print("j=", j)
        for i, rhoB in enumerate(rhoB_arr):
            rhoA = rho0 - rhoB

            if skip_neg_mu and eta >= get_eta_c(rhoA, rhoB, alpha=1):
                sigma_arr[j, i] = -1
                theta_arr[j, i] = -1
                continue
            M_gg_0, M_hh_0, f_bar = get_M_0(eta, rhoA, rhoB, alpha, K)
            M_gg_base, M_gh_base, M_hg_base, M_hh_base = get_M_cross_base(K)
            theta_max = None
            sigma_max = None
            for i_theta, theta in enumerate(theta_1D):
                qx = np.cos(theta) * q
                qy = np.sin(theta) * q
                M = assemble_matrix(M_gg_0, M_hh_0, M_gg_base, M_gh_base, M_hg_base, M_hh_base, qx, qy)
                sigma_re = np.max(linalg.eigvals(M).real)
                if theta_max is None or sigma_re > sigma_max:
                    theta_max = theta
                    sigma_max = sigma_re
            sigma_arr[j, i] = sigma_max
            theta_arr[j, i] = theta_max
    
    f_out = f"data/kinetic/r0{rho0:g}_q{q:.4f}_K{K:d}_n{n:d}_nt{n_theta}.npz"
    np.savez_compressed(f_out, rhoB=rhoB_arr, eta=eta_arr, sigma=sigma_arr, theta=theta_arr)


def cal_sigma_q_theta_fixed_rhoB(K=50, n=50, n_theta=30, skip_neg_mu=False):
    rhoA_arr = np.linspace(0, 1.2, n)
    eta_arr = np.linspace(0, 0.5, n)
    q = 1e-3
    rhoB = 0.03
    alpha = 1
    sigma_arr = np.zeros((eta_arr.size, rhoA_arr.size))
    theta_arr = np.zeros_like(sigma_arr)
    theta_1D = np.linspace(0, np.pi/2, n_theta)

    f_bar_left = None
    for i, rhoA in enumerate(rhoA_arr):
        print("i=", i)
        rho0 = rhoA + rhoB
        eta_c = get_eta_c(rhoA, rhoB, alpha=alpha)
        f_bar_down = None

        for j, eta in enumerate(eta_arr):
            if skip_neg_mu and eta >= eta_c:
                sigma_arr[j:, i] = -1
                theta_arr[j:, i] = -1
                break

            M_gg_0, M_hh_0, f_bar_down = get_M_0(eta, rhoA, rhoB, alpha, K, f_bar_pre=f_bar_down)

            M_gg_base, M_gh_base, M_hg_base, M_hh_base = get_M_cross_base(K)
            theta_max = None
            sigma_max = None
            for i_theta, theta in enumerate(theta_1D):
                qx = np.cos(theta) * q
                qy = np.sin(theta) * q
                M = assemble_matrix(M_gg_0, M_hh_0, M_gg_base, M_gh_base, M_hg_base, M_hh_base, qx, qy)

                sigma_re = np.max(linalg.eigvals(M).real)
                if theta_max is None or sigma_re > sigma_max:
                    theta_max = theta
                    sigma_max = sigma_re
            sigma_arr[j, i] = sigma_max
            theta_arr[j, i] = theta_max
    
    f_out = f"data/kinetic/rB{rhoB:g}_q{q:.4f}_K{K:d}_n{n:d}_nt{n_theta}.npz"
    np.savez_compressed(f_out, rhoA=rhoA_arr, eta=eta_arr, sigma=sigma_arr, theta=theta_arr)


def lin_diagram_w_color(fin, ax=None, ret_im=False):
    if ax is None:
        show_fig = True
        fig, ax = plt.subplots(1, 1, constrained_layout=True)
    else:
        show_fig = False

    with np.load(fin, "r") as data:
        if "rhoB" in data:
            x = data["rhoB"]
        else:
            x = data["rhoA"]
        eta = data["eta"]
        sigma = data["sigma"]
        state = np.zeros_like(sigma)
        mask = sigma > 1e-9
        state[mask] = 1
        theta = data["theta"]

    dx = x[1] - x[0]
    dy = eta[1] - eta[0]
    extent = [x[0]-dx/2, x[-1]+dx/2, eta[0]-dy/2, eta[-1]+dy/2]


    palette = plt.cm.plasma.with_extremes(over='w', under="tab:grey", bad='b')
    mask = np.logical_and(state == 0, theta >= 0)
    theta[mask] = np.pi

    im = ax.imshow(theta, origin="lower", cmap=palette, vmin=0, vmax=np.pi/2, extent=extent, aspect="auto")

    if show_fig:
        plt.show()
        plt.close()
    elif ret_im:
        return im


if __name__ == "__main__":
    K = 200

    # cal_sigma_q_theta_fixed_rhoB(K=2, n=200, n_theta=60, skip_neg_mu=True)
    cal_sigma_q_theta_fixed_rho0(K=1, n=200, n_theta=60, skip_neg_mu=True)



    # fins = ["data/kinetic/r01_q0.0010_K100_n200_nt60.npz", "data/kinetic/rB0.03_q0.0010_K100_n200_nt60.npz"]

    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), constrained_layout=True)
    # lin_diagram_w_color(fins[0], ax1)
    # lin_diagram_w_color(fins[1], ax2)
    # plt.show()
    # plt.close()
