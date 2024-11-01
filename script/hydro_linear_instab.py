import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy import linalg
from kinetic_linear_instab import get_P_k, get_eta_c, lin_diagram_w_color
import sys


def get_coeff(rA, rB, eta, alpha=1):
    P1 = get_P_k(1, eta)
    P2 = get_P_k(2, eta)
    rho0 = rA + rB
    mu = (2*alpha * (2*rA + rB)/ (np.pi * rho0) + 1) * P1 - (1+alpha)
    nu = 0.25 / (alpha+1-P2)
    gamma = 4 * alpha * nu / rho0 * (P2 - 2 * P1 / (3 * np.pi))
    kappa = 4 * alpha * nu / rho0 * (P2 + 2 * P1 / (3 * np.pi))
    xi = (4 * alpha / rho0)**2 * nu/(np.pi * 3) * P1 * P2
    mu_prime = 2 * alpha * P1 / np.pi * rB / rho0 **2
    xi_prime = (4 * alpha) ** 2 * nu / (np.pi * 3) * P1 * P2 * (-2) / rho0**3
    return mu, nu, gamma, kappa, xi, mu_prime, xi_prime


def get_M(q, phi, mu, nu, gamma, kappa, xi, mu_prime, xi_prime, N=3):
    M = np.zeros((N, N), complex)

    sqrt_mu_xi = np.sqrt(mu/xi)
    qx = q * np.cos(phi)
    qy = q * np.sin(phi)
    M[0, 1] = -1j * qx
    M[1, 0] = -0.5j * qx + (mu_prime -xi_prime * mu / xi) * sqrt_mu_xi
    M[1, 1] = -nu * q**2 - 1j * qx * sqrt_mu_xi * gamma - 2 * mu
    if N == 3:
        M[0, 2] = -1j * qy
        M[1, 2] = -1j * qy * sqrt_mu_xi * kappa
        M[2, 0] = -0.5j * qy
        M[2, 1] = 1j * qy * sqrt_mu_xi * kappa
        M[2, 2] = -nu * q**2 - 1j * qx * sqrt_mu_xi * gamma
    return M

def get_sigma(q, mu, nu, gamma, kappa, xi, mu_prime, xi_prime):
    sqrt_mu_xi = np.sqrt(mu/xi)
    a00 = 0
    a01 = -1j * q
    a10 = -0.5j * q + (mu_prime -xi_prime * mu / xi) * sqrt_mu_xi
    a11 = -nu * q**2 - 1j * q * sqrt_mu_xi * gamma - 2 * mu
    b = -(a00 + a11)
    c = a00 * a11 - a01 * a10
    Delta = b**2 - 4 * c
    root = 0.5 * (-b + np.sqrt(Delta))
    return root.real


def cal_sigma_q_theta_fixed_rho0(n=50, n_theta=30, skip_neg_mu=True, N=3):
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
            mu, nu, gamma, kappa, xi, mu_prime, xi_prime = get_coeff(rhoA, rhoB, eta, alpha=1)
            theta_max = None
            sigma_max = None
            for i_theta, theta in enumerate(theta_1D):
                M = get_M(q, theta, mu, nu, gamma, kappa, xi, mu_prime, xi_prime, N)
                try:
                    sigma_re = np.max(linalg.eigvals(M).real)
                except:
                    print(mu, nu, gamma, kappa, xi, mu_prime, xi_prime)
                    print(M)
                    sys.exit(1)
                if theta_max is None or sigma_re > sigma_max:
                    theta_max = theta
                    sigma_max = sigma_re
            sigma_arr[j, i] = sigma_max
            theta_arr[j, i] = theta_max

    if N == 2:
        f_out = f"data/hydros/r{rho0:g}_q{q:.4f}_n{n:d}_nt{n_theta}_N{N}.npz"
    else:
        f_out = f"data/hydros/r{rho0:g}_q{q:.4f}_n{n:d}_nt{n_theta}.npz"
    np.savez_compressed(f_out, rhoB=rhoB_arr, eta=eta_arr, sigma=sigma_arr, theta=theta_arr)


def cal_sigma_q_theta_fixed_rhoB(n=50, n_theta=30, skip_neg_mu=True):
    rhoA_arr = np.linspace(0, 1.2, n)
    eta_arr = np.linspace(0, 0.5, n)
    q = 1e-3
    rhoB = 0.03
    alpha = 1
    sigma_arr = np.zeros((eta_arr.size, rhoA_arr.size))
    theta_arr = np.zeros_like(sigma_arr)
    theta_1D = np.linspace(0, np.pi/2, n_theta)

    for i, rhoA in enumerate(rhoA_arr):
        print("i=", i)
        rho0 = rhoA + rhoB
        for j, eta in enumerate(eta_arr):

            mu, nu, gamma, kappa, xi, mu_prime, xi_prime = get_coeff(rhoA, rhoB, eta, alpha=1)
            theta_max = None
            sigma_max = None
            if skip_neg_mu and mu <= 0:
                sigma_arr[j:, i] = -1
                theta_arr[j:, i] = -1
                break
            theta_max = None
            sigma_max = None
            for i_theta, theta in enumerate(theta_1D):
                M = get_M(q, theta, mu, nu, gamma, kappa, xi, mu_prime, xi_prime)
                sigma_re = np.max(linalg.eigvals(M).real)
                if theta_max is None or sigma_re > sigma_max:
                    theta_max = theta
                    sigma_max = sigma_re
            sigma_arr[j, i] = sigma_max
            theta_arr[j, i] = theta_max
    
    f_out = f"data/hydros//rB{rhoB:g}_q{q:.4f}_n{n:d}_nt{n_theta}.npz"
    np.savez_compressed(f_out, rhoA=rhoA_arr, eta=eta_arr, sigma=sigma_arr, theta=theta_arr)


def func(rhoB, eta, rho0, alpha):
    mu, nu, gamma, kappa, xi, mu_prime, xi_prime = get_coeff(rho0-rhoB, rhoB, eta, alpha)
    left = (mu_prime - xi_prime * mu/xi - gamma * mu)**2
    right = gamma**2 * mu**2 + 2 * mu * xi
    return left - right


def find_eta_para():
    x_arr = np.linspace(0, 0.5, 400)
    y = 0.
    f = np.array([func(x, y, 1, 1) for x in x_arr])

    eta_arr = np.linspace(0, 0.5, 100)
    rhoB_arr = np.zeros_like(eta_arr)
    x0 = 0.3
    for i, eta in enumerate(eta_arr):
        sol = root_scalar(func, args=(eta, 1, 1), x0=x0, x1=0,  method="secant")
        x0 = sol.root
        rhoB_arr[i] = x0
    return rhoB_arr, eta_arr

def find_eat_thresh():
    rhoB_arr = np.linspace(0, 0.5, 400)
    rho0 = 1
    eta_c = get_eta_c(rho0-rhoB_arr, rhoB_arr, alpha=1)
    mask = eta_c > 0
    n = np.sum(mask) + 1
    return rhoB_arr[:n], eta_c[:n]


if __name__ == "__main__":
    # cal_sigma_q_theta_fixed_rho0(n=400, n_theta=60, skip_neg_mu=True, N=2)
    # cal_sigma_q_theta_fixed_rhoB(n=400, n_theta=60, skip_neg_mu=True)

    fins = ["data/hydros/r1_q0.0010_n400_nt60.npz",
            "data/hydros/r1_q0.0010_n400_nt60_N2.npz"
            # "data/hydros/rB0.03_q0.0010_n400_nt60.npz"
            ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), constrained_layout=True)

    lin_diagram_w_color(fins[0], ax1)
    lin_diagram_w_color(fins[1], ax2)

    ax1.set_ylim(0, 0.5)
    ax1.set_xlim(0, 0.5)

    plt.show()
    plt.close()


