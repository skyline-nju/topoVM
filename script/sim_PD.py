import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

def read_states(fin):
    with open(fin, "r") as f:
        lines = f.readlines()
        n = len(lines)
        x, y = np.zeros((2, n))
        state = np.zeros(n, int)
        for i, line in enumerate(lines):
            s = line.rstrip("\n").split()
            x[i] = float(s[0])
            y[i] = float(s[1])
            state[i] = int(s[2])
    return x, y, state


def rhoB_vs_eta(ax, ms=4):
    fin = "data/simulation_phase_diagram/rho0=1.dat"
    rhoB, eta, state = read_states(fin)

    mask = state == 0
    ax.plot(rhoB[mask], eta[mask], "o", c="tab:grey", ms=ms, fillstyle="none", label="disordered\ngas")

    mask = state == 1
    ax.plot(rhoB[mask], eta[mask], "x", c="tab:orange", ms=ms, fillstyle="none", label="polar liquid")

    mask = state == 2
    ax.plot(rhoB[mask], eta[mask], "s", c="tab:blue", ms=ms, fillstyle="none", label="bands")

    ax.set_xlim(0, 0.505)
    ax.set_ylim(0.045, 0.505)

    x1 = [0, 0.006, 0.008, 0.015, 0.023, 0.031, 0.038, 0.042, 0.04, 0.03, 0.02, 0.01, 0]
    y1 = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.37, 0.415, 0.442, 0.455, 0.465]
    ax.plot(x1, y1, c="k", lw=0.5)
    ax.fill_betweenx(y1, x1, 0, color="tab:orange", alpha=0.25)

    x2 = np.array([0.478, 0.475, 0.468, 0.45, 0.42, 0.38, 0.33, 0.27, 0.235, 0.12, 0.02])
    y2 = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.37, 0.415, 0.442]
    ax.plot(x2, y2, c="k", lw=0.5)
    ax.fill_betweenx(y2, x1[:x2.size], x2, color="tab:blue", alpha=0.25)

    
    y3 = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.37, 0.415, 0.442, 0.455, 0.465, 0.505]
    x3 = [0.478, 0.475, 0.468, 0.45, 0.42, 0.38, 0.33, 0.27, 0.235, 0.12, 0.02, 0.01, 0, 0]
    ax.fill_betweenx(y3, x3, 1, color="tab:grey", alpha=0.25)


def rhoA_vs_eta(ax, ms=4):
    fin = "data/simulation_phase_diagram/rhoB=0.03.dat"
    rhoA, eta, state = read_states(fin)

    mask = state == 0
    ax.plot(rhoA[mask], eta[mask], "o", c="tab:grey", ms=ms, fillstyle="none", label="disordered\ngas")

    mask = state == 1
    ax.plot(rhoA[mask], eta[mask], "x", c="tab:orange", ms=ms, fillstyle="none", label="polar liquid")

    mask = state == 2
    ax.plot(rhoA[mask], eta[mask], "s", c="tab:blue", ms=ms, fillstyle="none", label="bands")

    ax.set_xlim(0, 1.205)
    ax.set_ylim(0.045, 0.505)

    x1 = np.array([0.004, 0.008, 0.02, 0.025, 0.031, 0.038, 0.042, 0.04, 0.03, 0.02, 0.01])
    y1 = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.37, 0.415, 0.442, 0.455]
    x1 = 0.03/x1-0.03
    # print(x1)
    ax.plot(x1, y1,"-", c="k", lw=0.5)
    ax.fill_betweenx(y1, x1, 2, color="tab:orange", alpha=0.25)

    x2 = np.array([0.33, 0.32, 0.3, 0.28,
                   0.26, 0.24, 0.2, 0.12,
                   0.10, 0.06, 0.02])
    y2 = np.array([0, 0.05, 0.1, 0.15,
                   0.2, 0.25, 0.3, 0.35,
                   0.37, 0.415, 0.442])

    x2 = 0.03/x2 - 0.03
    ax.plot(x2, y2,"-", c="k", lw=0.5)
    x1_new = np.zeros(x1.size + 1)
    x1_new[0] = 0.00001
    x1_new[1:] = x1
    ax.fill_betweenx(y2, x2, x1_new[:x2.size], color="tab:blue", alpha=0.25)

    y2_new = np.zeros(y2.size +1)
    y2_new[:-1] = y2
    y2_new[-1] = 1
    x2_new = np.zeros(y2.size + 1)
    x2_new[:-1] = x2
    x2_new[-1] = 10
    ax.fill_betweenx(y2_new, x2_new, 0, color="tab:grey", alpha=0.25)

if __name__ == "__main__":



    fig, (ax1, ax2) = plt.subplots(1, 2, constrained_layout=True, figsize=(8, 4), sharey=True)

    rhoB_vs_eta(ax1, ms=4)
    rhoA_vs_eta(ax2, ms=4)

    plt.show()
    plt.close()