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


def plot_band_profiles_varied_rhoB(axes=None, rhoB_arr=None, lw=1):
    if axes is None:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        flag_show = True
    else:
        ax1, ax2 = axes
        flag_show = False

    prefix = f"{root_sohrab}/topoVM/dissenters/L400_new2"
    if flag_show:
        coarse_grain_all(prefix, dx=4, pat="L400_400_d*_e0.100_r1_s*.gsd")

    npz_folder = f"{prefix}/coarse_grain_dx4"

    if rhoB_arr is None:
        rhoB_arr = [0.02, 0.03, 0.05, 0.1, 0.2, 0.3]
    seed = {0.02: 2097, 0.05: 2076, 0.2: 2044, 0.03:2073, 0.1:2037, 0.3:2049}
    
    for i, rhoB in enumerate(rhoB_arr):
        basename = f"L400_400_d{rhoB:.4f}_e0.100_r1_s{seed[rhoB]:d}.npz"
        fout = f"data/time_ave_profile/L400_400_rB{rhoB:.4f}_e0.1_r1.npz"
        fname = f"{npz_folder}/{basename}"
        if flag_show:
            x, rho_m, mx_m = get_mean_profiles(fname, 0, rho_threshold=4.5)
        else:
            with np.load(fout, "r") as data:
                x, rho_m, mx_m = data["x"], data["rhoA"], data["mA"]
        ax1.plot(x, rho_m, lw=lw, label=r"$%g$" % rhoB)
        ax2.plot(x, mx_m/rho_m, lw=lw)

        if flag_show:
            np.savez_compressed(fout, x=x, rhoA=rho_m, mA=mx_m)
    if flag_show:
        plt.show()
        plt.close()


def plot_band_profiles_varied_eta(axes=None, lw=1):
    if axes is None:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        flag_show = True
    else:
        ax1, ax2 = axes
        flag_show = False
    prefix = f"{root_sohrab}/topoVM/dissenters/L400_new2"
    if flag_show:
        coarse_grain_all(prefix, dx=4, pat="L400_400_d0.3000_e*_r1_s*.gsd")


    npz_folder = f"{prefix}/coarse_grain_dx4"

    eta_arr = [0.1, 0.2, 0.3]
    c_arr = ["tab:red", "tab:purple", "tab:cyan"]
    seed = 2049
    
    for i, eta in enumerate(eta_arr):
        basename = f"L400_400_d0.3000_e{eta:.3f}_r1_s{seed:d}.npz"
        fout = f"data/time_ave_profile/L400_400_rB0.3000_e{eta:.2f}_r1.npz"
        fname = f"{npz_folder}/{basename}"
        if flag_show:
            x, rho_m, mx_m = get_mean_profiles(fname, 0, rho_threshold=4)
        else:
            with np.load(fout, "r") as data:
                x, rho_m, mx_m = data["x"], data["rhoA"], data["mA"]
        ax1.plot(x, rho_m, lw=lw, label=r"$%g$" % eta, c=c_arr[i])
        ax2.plot(x, mx_m/rho_m, lw=lw, c=c_arr[i])

        if flag_show:
            np.savez_compressed(fout, x=x, rhoA=rho_m, mA=mx_m)
    if flag_show:
        plt.show()
        plt.close()


def get_instant_profiles():
    prefix = f"{root_sohrab}//topoVM/dissenters/L400_new"
    npz_folder = f"{prefix}/coarse_grain_dx4"
    basenames = [
        "L400_400_d0.0300_e0.100_r1_s3100.npz",
    ]

    iframes = [80]
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
    # plot_band_profiles_varied_rhoA()

    # get_instant_profiles()
    # coarse_grain_all(prefix, dx=4)

    plot_band_profiles_varied_rhoB()



    
        


