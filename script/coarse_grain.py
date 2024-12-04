import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from gsd import fl

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_rudabeh = "/run/user/1148/gvfs/sftp:host=rudabeh002,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def coarse_grain_one_frame(i_frame: int,
                           fgsd: fl.GSDFile,
                           fields: np.ndarray,
                           dx: float,
                           Lx: float,
                           Ly: float):
    ncols = int(Lx / dx)
    nrows = int(Ly / dx)
    bin_area = dx ** 2

    par_num = np.zeros((2, nrows, ncols), int)
    mx, my = np.zeros((2, 2, nrows, ncols))

    pos = fgsd.read_chunk(frame=i_frame, name="particles/position")
    try:
        id_type = fgsd.read_chunk(frame=i_frame, name="particles/typeid")
    except KeyError:
        id_type = np.zeros(pos.shape[0], np.int32)

    x = pos[:,0] + Lx / 2
    y = pos[:,1] + Ly / 2
    ux = np.cos(pos[: ,2])
    uy = np.sin(pos[:, 2])

    col = (x/dx).astype(int)
    row = (y/dx).astype(int)
    col[col < 0] += ncols
    col[col >= ncols] -= ncols
    row[row < 0] += nrows
    row[row >= nrows] -= nrows
    
    for j in range(x.size):
        my_id = id_type[j]
        my_row = row[j]
        my_col = col[j]
        par_num[my_id, my_row, my_col] += 1
        mx[my_id, my_row, my_col] += ux[j]
        my[my_id, my_row, my_col] += uy[j]
    fields[i_frame][0:2] = par_num / bin_area
    fields[i_frame][2:4] = mx / bin_area
    fields[i_frame][4:6] = my / bin_area


def coarse_grain(fin: str, fout: str, dx: float):
    if os.path.exists(fout):
        with np.load(fout, "r") as data:
            t0 = data["t"]
            x0 = data["x"]
            y0 = data["y"]
            fields0 = data["fields"]
            fields_name0 = data["fields_name"]
            existed_frames = t0.size
    else:
        existed_frames = 0
    with fl.open(name=fin, mode="rb") as fgsd:
        if fgsd.nframes == 0:
            print("Warning, zero frame found in", fin)
        else:
            nframes = fgsd.nframes
            print(nframes, "frames found for", fin)
            Lx, Ly = fgsd.read_chunk(frame=0, name="configuration/box")[:2]
            ncols = int(Lx / dx)
            nrows = int(Ly / dx)
            x = np.arange(ncols) * dx + dx / 2
            y = np.arange(nrows) * dx + dx / 2
            
            t = np.zeros(nframes, int)

            fields = np.zeros((nframes, 6, nrows, ncols), np.single)

            if existed_frames > 0:
                t[:existed_frames] = t0
                fields[:existed_frames] = fields0

            for i_frame in range(existed_frames, nframes):
                try:
                    t[i_frame] = fgsd.read_chunk(frame=i_frame, name="configuration/step")[0]
                except KeyError:
                    t[i_frame] = 0
                    print("Warning, not found t[%d]" % i_frame)
                print("coarse-graining frame %d/%d, t=%g..." % (i_frame, nframes, t[i_frame]))
                coarse_grain_one_frame(i_frame, fgsd, fields, dx, Lx, Ly)
            
            fields_name = np.array(["rho_A", "rho_B", "mx_A", "mx_B", "my_A", "my_B"])

            np.savez_compressed(fout, t=t, x=x, y=y, fields=fields, fields_name=fields_name)


def coarse_grain_all(prefix, dx, forced_updating=False, pat="*.gsd"):
    out_folder = f"{prefix}/coarse_grain_dx{dx:g}"
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    
    fins = glob.glob(f"{prefix}/{pat}")
    for fin in fins:
        basename = os.path.basename(fin)
        out_basename = basename.replace(".gsd", ".npz")
        fout = f"{out_folder}/{out_basename}"
        if os.path.exists(fout):
            if not forced_updating and os.path.getmtime(fout) > os.path.getmtime(fin):
                print("skipping", fin)
                continue
        coarse_grain(fin, fout, dx)


if __name__ == "__main__":
    from datetime import datetime

    prefix = f"{root_sohrab}/topoVM/dissenters/L2000"
    coarse_grain_all(prefix, dx=4, pat="L2000_2000_d0.0000_e0.440_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=2, pat="L2000_2000_d0.0000_e0.440_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=40, pat="L2000_2000_d0.0000_e0.450_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=4, pat="L2000_2000_d0.0300_e0.100_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=4, pat="L2000_2000_d0.0300_e0.300_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=2, pat="L2000_2000_d0.0300_e0.300_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=4, pat="L2000_2000_d0.0300_e0.420_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=4, pat="L2000_2000_d0.0000_e0.480_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=4, pat="L2000_2000_d0.0000_e0.500_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=25, pat="L2000_2000_d0.0300_e0.100_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=25, pat="L2000_2000_d0.0300_e0.300_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=25, pat="L2000_2000_d0.0300_e0.420_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=10, pat="L2000_2000_d0.0300_e0.300_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=10, pat="L2000_2000_d0.0300_e0.420_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=20, pat="L2000_2000_d0.0300_e0.300_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=20, pat="L2000_2000_d0.0300_e0.420_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=40, pat="L2000_2000_d0.0300_e0.300_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=40, pat="L2000_2000_d0.0300_e0.420_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=40, pat="L2000_2000_d0.0000_e0.440_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=80, pat="L2000_2000_d0.0300_e0.300_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=80, pat="L2000_2000_d0.0300_e0.420_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=100, pat="L2000_2000_d0.0300_e0.300_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=100, pat="L2000_2000_d0.0300_e0.420_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=125, pat="L2000_2000_d0.0300_e0.300_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=125, pat="L2000_2000_d0.0300_e0.420_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=200, pat="L2000_2000_d0.0300_e0.300_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=200, pat="L2000_2000_d0.0300_e0.420_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=400, pat="L2000_2000_d0.0300_e0.300_r1_s3100.gsd")

    prefix = f"{root_sohrab}/topoVM/dissenters/L1000"
    coarse_grain_all(prefix, dx=4, pat="L1000_1000_d0.0300_e0.300_r1_s3100.gsd")
    coarse_grain_all(prefix, dx=4, pat="L1000_1000_d0.0300_e0.450_r1_s3100.gsd")

    prefix = f"{root_sohrab}/topoVM/dissenters/L400_new"
    coarse_grain_all(prefix, dx=40, pat="*.gsd")
    coarse_grain_all(prefix, dx=4, pat="L400_400_d0.2000_e0.100_r1_s3100.gsd")

    prefix = f"{root_sohrab}/topoVM/dissenters/L200_new"
    coarse_grain_all(prefix, dx=40, pat="*.gsd")

    prefix = f"{root_sohrab}/topoVM/dissenters/L800_new2"
    coarse_grain_all(prefix, dx=40, pat="*.gsd")

    prefix = f"{root_sohrab}/topoVM/dissenters/L1600"
    coarse_grain_all(prefix, dx=40, pat="*.gsd")

    # prefix = f"{root_sohrab}/topoVM/dissenters/L2000"
    # coarse_grain_all(prefix, dx=200, pat="*_d0.0300_*.gsd")

    # prefix = f"{root_sohrab}/topoVM/dissenters/L2000"
    # coarse_grain_all(prefix, dx=125, pat="*_d0.0300_*.gsd")

    # prefix = f"{root_sohrab}/topoVM/dissenters/L800_varied_rhoA"
    # coarse_grain_all(prefix, dx=4, pat="*.gsd")