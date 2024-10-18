import numpy as np
import matplotlib.pyplot as plt
from coarse_grain import coarse_grain_all
root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"


def find_interface(rho_arr, rho_threshold=2):
    xs = []
    n = rho_arr.size
    for i in range(n):
        if rho_arr[i-1] > rho_threshold and rho_arr[i] <= rho_threshold:
            is_valid = True
            for j in range(20):
                i2 = i + j
                if i2 >= n:
                    i2 -= n
                if rho_arr[i2] > rho_threshold:
                    is_valid = False
                    break
            if is_valid:
                for j in range(1, 20):
                    i2 = i - j
                    print(i, i2)
                    if rho_arr[i2] < rho_threshold:
                        is_valid = False
                        break
            if is_valid:
                xs.append(i)
    if len(xs) == 0:
        for i in range(n):
            if rho_arr[i-1] > rho_threshold and rho_arr[i] <= rho_threshold:
                is_valid = True
                for j in range(20):
                    i2 = i + j
                    if i2 >= n:
                        i2 -= n
                    if rho_arr[i2] > rho_threshold:
                        is_valid = False
                        break
                if is_valid:
                    xs.append(i)
    return xs


def get_mean_profiles(fname, ncut, rho_threshold=2):
    with np.load(fname, "r") as data:
        t = data["t"]
        x = data["x"]
        fields = data["fields"]

        rhoA = np.mean(fields[:, 0], axis=1)
        # rhoB = fields[:, 1]
        mA_x = np.mean(fields[:, 2], axis=1)
        # mA_y = np.mean(fields[:, 4], axis=1)
        px = mA_x / rhoA
    
    rho_shift = np.zeros_like(rhoA)
    m_shift = np.zeros_like(mA_x)
    p_shift = np.zeros_like(px)
    for i in range(t.size):
        xs = find_interface(rhoA[i], rho_threshold=rho_threshold)
        if len(xs) > 1 or len(xs) == 0:
            fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, constrained_layout=True)
            ax1.plot(x, rhoA[i], "-o")
            ax2.plot(x, px[i])
            ax1.axhline(2, linestyle="dashed", c="tab:red")
            for idx_x in xs:
                ax1.axvline(x[idx_x], linestyle="dotted", c="tab:red")
                ax2.axvline(x[idx_x], linestyle="dotted", c="tab:red")
            plt.title(xs)
            plt.show()
            plt.close()
        elif len(xs) == 1:
            idx = xs[0]
            shift = 150 - idx
            rhoA_new = np.roll(rhoA[i], shift)
            rho_shift[i] = rhoA_new
            m_shift[i] = np.roll(mA_x[i], shift)

    rho_m = np.mean(rho_shift[ncut:], axis=0)
    mx_m = np.mean(m_shift[ncut:], axis=0)
    return x, rho_m, mx_m


def plot_band_profiles_varied_rhoA():
    prefix = f"{root_sohrab}//topoVM/dissenters/L800_varied_rhoA"
    coarse_grain_all(prefix, dx=4)

    npz_folder = f"{prefix}/coarse_grain_dx4"

    rho0_arr = [1, 1.2, 1.4, 1.6, 1.8]
    m_arr = []

    basenames = ["L800_800_d0.1000_e0.300_r1_s2011.npz",
                 "L800_800_d0.083333_e0.300_r1.2_s2001.npz",
                 "L800_800_d0.071429_e0.300_r1.4_s2001.npz",
                 "L800_800_d0.062500_e0.300_r1.6_s2001.npz",
                 "L800_800_d0.055556_e0.300_r1.8_s2001.npz"
                 ]
    ncuts = [60, 30, 30, 30, 30, 20]
    
    for i, basename in enumerate(basenames):
        fname = f"{npz_folder}/{basename}"
        x, rho_m, mx_m = get_mean_profiles(fname, 20)
        plt.plot(x, rho_m)
        m_arr.append(np.mean(mx_m))

        fout = f"data/time_ave_profile/L800_800_rB0.1_e0.3_r{rho0_arr[i]:.2f}.npz"
        np.savez_compressed(fout, x=x, rhoA=rho_m, mA=mx_m)
    plt.show()
    plt.close()

    plt.plot(rho0_arr, m_arr, "-o")
    plt.show()
    plt.close()


def get_instant_profiles():
    prefix = f"{root_sohrab}//topoVM/dissenters/L800_varied_rhoA"
    npz_folder = f"{prefix}/coarse_grain_dx4"
    basenames = [
        "L800_800_d0.1000_e0.300_r1_s2004.npz",
        "L800_800_d0.1000_e0.300_r1_s2002.npz",
        "L800_800_d0.1000_e0.300_r1_s2011.npz",
    ]

    iframes = [-1, -1, 197]
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, constrained_layout=True)
    for i, basename in enumerate(basenames):
        fname = f"{npz_folder}/{basename}"
        with np.load(fname, "r") as data:
            t = data["t"]
            x = data["x"]
            fields = data["fields"]

            j = iframes[i]
            rhoA = np.mean(fields[j, 0], axis=0)
            # rhoB = fields[:, 1]
            mA_x = np.mean(fields[j, 2], axis=0)
            # mA_y = np.mean(fields[:, 4], axis=)
            px = mA_x / rhoA
            ax1.plot(x, rhoA)
            ax2.plot(x, px)
        fout = f"data/profile/{basename}"
        np.savez_compressed(fout, x=x, rhoA=rhoA, mA=mA_x)
    plt.show()
    plt.close()

        

if __name__ == "__main__":
    plot_band_profiles_varied_rhoA()

    # get_instant_profiles()
    # prefix = f"{root_sohrab}/topoVM/dissenters/L400_new"
    # coarse_grain_all(prefix, dx=4)

    # npz_folder = f"{prefix}/coarse_grain_dx10"
    
    # basename = "L400_400_d0.0300_e0.100_r1_s3100.npz"
    # fname = f"{npz_folder}/{basename}"
    # with np.load(fname, "r") as data:
    #     t = data["t"]
    #     x = data["x"]
    #     fields = data["fields"]

    #     rhoA = np.mean(fields[:, 0], axis=1)
    #     # rhoB = fields[:, 1]
    #     mA_x = np.mean(fields[:, 2], axis=1)
    #     # mA_y = np.mean(fields[:, 4], axis=1)
    #     px = mA_x / rhoA
    
    # for i in range(t.size):
    #     fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, constrained_layout=True)
    #     ax1.plot(x, rhoA[i])
    #     ax2.plot(x, px[i])
    #     plt.suptitle(r"$t=%g$" % t[i])
    #     plt.show()
    #     plt.close()

    
        


