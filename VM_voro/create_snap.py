import sys
import numpy as np
from gsd import hoomd
from read_gsd import read_one_frame


def sort_by_type(snap) -> hoomd.Snapshot:
    N = snap.particles.N
    pos = np.zeros((N, 3), dtype=np.float32)
    type_id = np.zeros(N, dtype=np.uint32)

    mask_aligner = snap.particles.typeid == 0
    mask_dissenter = snap.particles.typeid == 1
    n_dis = np.sum(mask_dissenter)
    type_id[:n_dis] = 1
    pos[:n_dis] = snap.particles.position[mask_dissenter, :]
    pos[n_dis:] = snap.particles.position[mask_aligner, :]
    snap.particles.position = pos
    snap.particles.typeid = type_id
    return snap
    

def duplicate(s: hoomd.Snapshot, nx: int, ny: int) -> hoomd.Snapshot:
    N = s.particles.N * nx * ny
    lx = s.configuration.box[0]
    ly = s.configuration.box[1]
    Lx, Ly = lx * nx, ly * ny
    pos = np.zeros((N, 3), dtype=np.float32)
    type_id = np.zeros(N, dtype=np.uint32)
    for j in range(ny):
        for i in range(nx):
            beg = (j * nx + i) * s.particles.N
            end = beg + s.particles.N
            pos[beg:end, 0] = s.particles.position[:, 0] + lx / 2 + i * lx
            pos[beg:end, 1] = s.particles.position[:, 1] + ly / 2 + j * ly
            pos[beg:end, 2] = s.particles.position[:, 2]
            type_id[beg:end] = s.particles.typeid
    pos[:, 0] -= Lx / 2
    pos[:, 1] -= Ly / 2
    s2 = hoomd.Snapshot()
    s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s2.particles.N = N
    s2.particles.position = pos
    s2.particles.typeid = type_id
    s2.particles.types = s.particles.types
    s2.configuration.step = 0
    return s2


def scale(s: hoomd.Snapshot, nx: int, ny: int, eps=0) -> hoomd.Snapshot:
    lx = s.configuration.box[0]
    ly = s.configuration.box[1]
    Lx, Ly = lx * nx, ly * ny
    if isinstance(nx, int) and isinstance(ny, int):
        N = s.particles.N * nx * ny
        pos = np.zeros((N, 3), dtype=np.float32)
        type_id = np.zeros(N, dtype=np.uint32)
        for i in range(nx * ny):
            beg = i * s.particles.N
            end = beg + s.particles.N
            pos[beg:end, 0] = s.particles.position[:, 0] * nx
            pos[beg:end, 1] = s.particles.position[:, 1] * ny
            pos[beg:end, 2] = s.particles.position[:, 2]
            type_id[beg:end] = s.particles.typeid
        if nx > 1:
            pos[:, 0] += (np.random.rand(N) - 0.5) * eps * nx
            mask = pos[:, 0] < Lx/2
            pos[:, 0][mask] += Lx
            mask = pos[:, 0] >= Lx/2
            pos[:, 0][mask] -= Lx
        if ny > 1:
            pos[:, 1] += (np.random.rand(N) - 0.5) * eps * ny
            mask = pos[:, 1] < Ly/2
            pos[:, 1][mask] += Ly
            mask = pos[:, 1] >= Ly/2
            pos[:, 1][mask] -= Ly
        s2 = hoomd.Snapshot()
        s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
        s2.particles.N = N
        s2.particles.position = pos
        s2.particles.typeid = type_id
        s2.particles.types = s.particles.types
        s2.configuration.step = 0
    else:
        s.particles.position[:, 0] *= nx
        s.particles.position[:, 1] *= ny
        s.configuration.box = [Lx, Ly, 1, 0, 0, 0]

        mask_A = s.particles.typeid == 0
        mask_B = s.particles.typeid == 1
        rho_A = np.sum(mask_A) / (lx * ly)
        rho_B = np.sum(mask_B) / (lx * ly)
        print("rho_A=", rho_A, "rho_B=", rho_B)
        s2 = adjust_density(s, rho_A, rho_B, mode="copy")
    return s2


def adjust_density(s: hoomd.Snapshot, phi_A: float, phi_B: float, mode="copy", R=5) -> hoomd.Snapshot:
    def add_new_particles(n):
        pos_new = np.zeros((n, 3), dtype=np.float32)
        pos_new[:, 0] = (np.random.rand(n) - 0.5) * Lx
        pos_new[:, 1] = (np.random.rand(n) - 0.5) * Ly
        pos_new[:, 2] = (np.random.rand(n) - 0.5) * np.pi * 2
        return pos_new


    def add_new_particles_locally(n, R, shape="circle", Ly=None):
        pos_new = np.zeros((n, 3), dtype=np.float32)
        # r = np.random.rand(n) * R
        # theta = np.random.rand(n) * np.pi * 2
        # pos_new[:, 0] = r * np.cos(theta)
        # pos_new[:, 1] = r * np.sin(theta)
        if shape == "circle":
            i = 0
            while i < n:
                x, y = np.random.rand(2) * 2 * R - R
                if x**2 + y**2 < R**2:
                    pos_new[i, 0] = x
                    pos_new[i, 1] = y
                    i += 1
        elif shape == "rect":
            if Ly is None:
                print("Error, Ly is needed when using rect shape")
                sys.exit()
            pos_new[:, 0] = np.random.rand(n) * 2 * R - R
            pos_new[:, 1] = np.random.rand(n) * Ly - Ly / 2
        pos_new[:, 2] = (np.random.rand(n) - 0.5) * np.pi * 2
        return pos_new

    def extend_pos(pos_new, pos_old, mode):
        n_new = pos_new.shape[0]
        n_old = pos_old.shape[0]
        m = n_new // n_old
        n_res = n_new - m * n_old
        end = 0
    
        for i in range(m):
            beg = i * n_old
            end = beg + n_old
            pos_new[beg: end] = pos_old
        
        if n_res > 0:
            if mode == "rand":
                pos_new[end:] = add_new_particles(n_res)
            elif mode == "copy":
                np.random.shuffle(pos_old)
                pos_new[end:] = pos_old[:n_res]
            elif mode == "local":
                pos_new[end:] = add_new_particles_locally(n_res, R)
            elif mode == "local_rect":
                pos_new[end:] = add_new_particles_locally(n_res, R, shape="rect", Ly=Ly)
            else:
                print("Error, mode should be copy or rand")
        

    Lx = int(s.configuration.box[0])
    Ly = int(s.configuration.box[1])
    print(Lx, Ly)

    n_A = int(round(Lx * Ly * phi_A))
    n_B = int(round(Lx * Ly * phi_B))
    print("new:", n_A+n_B, n_A, n_B)

    pos_A = np.zeros((n_A, 3), dtype=np.float32)
    pos_B = np.zeros((n_B, 3), dtype=np.float32)

    mask_A = s.particles.typeid == 0
    n_A0 = np.sum(mask_A)
    mask_B = s.particles.typeid == 1
    n_B0 = np.sum(mask_B)
    print("old:", s.particles.N, n_A0, n_B0)

    pos_A0 = s.particles.position[mask_A]
    if n_A0 < n_A:
        extend_pos(pos_A, pos_A0, mode)
        print("add %d A particles" % (n_A - n_A0))
    elif n_A0 > n_A:
        np.random.shuffle(pos_A0)
        pos_A = pos_A0[:n_A]
        print("remove %d A particles" % (n_A0 - n_A))
    else:
        pos_A = pos_A0
    
    pos_B0 = s.particles.position[mask_B]
    if n_B0 < n_B:
        extend_pos(pos_B, pos_B0, mode)
        print("add %d B particles" % (n_B - n_B0))
    elif n_B0 > n_B:
        np.random.shuffle(pos_B0)
        pos_B = pos_B0[:n_B]
        print("remove %d B particles" % (n_B0 - n_B))
    else:
        pos_B = pos_B0
    
    s2 = hoomd.Snapshot()
    s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s2.particles.N = n_A + n_B
    s2.particles.position = np.vstack((pos_A, pos_B))
    s2.particles.typeid = np.zeros(s2.particles.N, dtype=np.float32)
    s2.particles.typeid[n_A:] = 1
    s2.particles.types = s.particles.types
    s2.configuration.step = 0
    return s2


def add_A_particles(s: hoomd.Snapshot, nA_added) -> hoomd.Snapshot:
    Lx = int(s.configuration.box[0])
    Ly = int(s.configuration.box[1])

    n_old = s.particles.N
    n_new = s.particles.N + nA_added
    s.particles.N = n_new
    pos = np.zeros((n_new, 3), dtype=np.float32)
    pos[:n_old] = s.particles.position
    type_id = np.zeros(n_new, dtype=np.uint32)
    type_id[:n_old] = s.particles.typeid

    pos[n_old:, 0] = (np.random.rand(nA_added) - 0.5) * Lx
    pos[n_old:, 1] = (np.random.rand(nA_added) - 0.5) * Ly
    pos[n_old:, 2] = (np.random.rand(nA_added) - 0.5) * np.pi * 2

    s2 = hoomd.Snapshot()
    s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s2.particles.N = n_new
    s2.particles.position = pos
    s2.particles.typeid = type_id
    s2.particles.types = s.particles.types
    s2.configuration.step = 0
    return s2




def inverse(x, y, theta, xc, yc):
    x_inv = 2 * xc - x
    y_inv = 2 * yc - y
    theta_inv = theta + np.pi
    theta_inv[theta_inv > np.pi] -= np.pi * 2
    theta_inv[theta_inv < -np.pi] += np.pi * 2
    return x_inv, y_inv, theta_inv


def flip_half(s):
    x, y, theta = s.particles.position.T
    Lx = s.configuration.box[0]
    mask1 = x < 0
    mask2 = x >= 0
    x1, y1, theta1 = x[mask1], y[mask1], theta[mask1]
    x2, y2, theta2 = x[mask2], y[mask2], theta[mask2]
    x1, y1, theta1 = inverse(x1, y1, theta1, -Lx / 4, 0)
    type_id1 = s.particles.typeid[mask1]
    type_id2 = s.particles.typeid[mask2]

    x = np.hstack((x1, x2))
    y = np.hstack((y1, y2))
    theta = np.hstack((theta1, theta2))
    type_id = np.hstack((type_id1, type_id2))
    s.particles.position = np.array([x, y, theta], dtype=np.float32).T
    s.particles.typeid = np.array(type_id, dtype=np.uint32)
    return s


def create_nematic_patten(s, nb=2):
    Lx = s.configuration.box[0]
    s = flip_half(s)
    if nb == 4:
        x = s.particles.position[:, 0]
        mask1 = np.logical_and(x >= -Lx/4, x < 0)
        mask2 = np.logical_and(x >= 0, x <= Lx/4)
        x[mask1] += Lx / 4
        x[mask2] -= Lx / 4
        s.particles.position[:, 0] = x
    s.configuration.step = 0
    return s


def create_flat_band(Lx, Ly, phi_A, phi_B, rho_A, rho_B, fname):
    s = hoomd.Snapshot()
    s.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s.configuration.step = 0
    nA = int(Lx * Ly * phi_A)
    nB = int(Lx * Ly * phi_B)
    n_tot = nA + nB
    s.particles.N = n_tot

    s.particles.types = ['A', 'B']
    s.particles.typeid = np.zeros(n_tot, np.uint32)
    s.particles.typeid[nA:] = 1

    pos = np.zeros((n_tot, 3), np.float32)
    l_A = Lx * (phi_A - rho_A[1]) / (rho_A[0] - rho_A[1])
    l_B = Lx * (phi_B - rho_B[1]) / (rho_B[0] - rho_B[1])
    nA_lt = int(l_A * Ly * rho_A[0])
    nB_lt = int(l_B * Ly * rho_B[0])
    nA_rt = nA - nA_lt
    nB_rt = nB - nB_lt
    pos[:nA_lt, 0] = np.random.rand(nA_lt) * l_A
    pos[nA_lt: nA, 0] = np.random.rand(nA_rt) * (Lx - l_A) + l_A
    pos[nA: nA+nB_lt, 0] = np.random.rand(nB_lt) * l_B
    pos[nA+nB_lt:, 0] = np.random.rand(nB_rt) * (Lx - l_B) + l_B 
    pos[:, 0] -= 0.5 * Lx
    pos[:, 1] = (np.random.rand(n_tot) - 0.5) * Ly
    pos[:, 2] = (np.random.rand(n_tot) - 0.5) * np.pi * 2
    s.particles.position = pos

    with hoomd.open(name=fname, mode="wb") as f:
        f.append(s)


def shift_v_slice(s: hoomd.Snapshot, x1: float, x2: float, dy: float) -> hoomd.Snapshot:
    Ly = s.configuration.box[1]
    x = s.particles.position[:, 0]
    y = s.particles.position[:, 1]

    mask = np.logical_and(x >= x1, x < x2)
    y[mask] += dy
    y[y >= Ly / 2] -= Ly
    y[y < -Ly/2] += Ly
    s.particles.position[:, 1] = y
    return s


def shift_h_slice(s: hoomd.Snapshot, y1: float, y2: float, dx: float) -> hoomd.Snapshot:
    Lx = s.configuration.box[0]
    x = s.particles.position[:, 0]
    y = s.particles.position[:, 1]

    mask = np.logical_and(y >= y1, y < y2)
    x[mask] += dx
    x[x >= Lx/2] -= Lx
    x[x < -Lx/2] += Lx
    s.particles.position[:, 0] = x
    return s

def shift_pos(s: hoomd.Snapshot, dx: float, dy: float) -> hoomd.Snapshot:
    Ly = s.configuration.box[1]
    Lx = s.configuration.box[0]

    x = s.particles.position[:, 0]
    y = s.particles.position[:, 1]

    x += dx
    y += dy
    if dy > 0:
        y[y >= Ly / 2] -= Ly
    else:
        y[y < -Ly / 2] += Ly
    if dx > 0:
        x[x >= Lx / 2] -= Lx
    else:
        x[x < -Lx / 2] += Lx
    s.particles.position[:, 0] = x
    s.particles.position[:, 1] = y
    return s


def shift_to_center(s: hoomd.Snapshot) -> hoomd.Snapshot:
    Ly = s.configuration.box[1]
    Lx = s.configuration.box[0]

    x = s.particles.position[:, 0] / Lx * np.pi * 2
    y = s.particles.position[:, 1] / Ly * np.pi * 2

    x_c = np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))) / (2 * np.pi) * Lx
    y_c = np.arctan2(np.mean(np.sin(y)), np.mean(np.cos(y))) / (2 * np.pi) * Ly

    return shift_pos(s, -x_c, -y_c)

def swap_h(s: hoomd.Snapshot, x1: float, x2: float, dx: float) -> hoomd.Snapshot:
    Lx = s.configuration.box[0]
    x = s.particles.position[:, 0]

    mask_a = np.logical_and(x >= x1, x < x2)
    mask_b = np.logical_and(x >= x1 + dx, x < x2 + dx)
    x[mask_a] += dx
    x[mask_b] -= dx
    x[x>= Lx/2] -= Lx/2
    x[x< -Lx/2] += Lx/2
    s.particles.position[:, 0] = x
    return s


def hstack(s1: hoomd.Snapshot, s2: hoomd.Snapshot) -> hoomd.Snapshot:
    lx1 = s1.configuration.box[0]
    ly1 = s1.configuration.box[1]
    lx2 = s2.configuration.box[0]
    ly2 = s2.configuration.box[1]
    if ly1 != ly2:
        print("Error, Ly=%g,%g" % (ly1, ly2))
        sys.exit(1)
    n1 = s1.particles.N
    N = n1 + s2.particles.N
    Lx, Ly = lx1 + lx2, ly1
    pos = np.zeros((N, 3), dtype=np.float32)
    type_id = np.zeros(N, dtype=np.uint32)

    dx1 = -lx2/2
    dx2 = lx1/2
    pos[:n1, 0] = s1.particles.position[:, 0] + dx1
    pos[n1:, 0] = s2.particles.position[:, 0] + dx2
    pos[:n1, 1] = s1.particles.position[:, 1]
    pos[n1:, 1] = s2.particles.position[:, 1]
    pos[:n1, 2] = s1.particles.position[:, 2]
    pos[n1:, 2] = s2.particles.position[:, 2]
    type_id[:n1] = s1.particles.typeid
    type_id[n1:] = s2.particles.typeid

    s2 = hoomd.Snapshot()
    s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s2.particles.N = N
    s2.particles.position = pos
    s2.particles.typeid = type_id
    s2.particles.types = s1.particles.types
    s2.configuration.step = 0
    return s2


def vstack(s1: hoomd.Snapshot, s2: hoomd.Snapshot) -> hoomd.Snapshot:
    lx1 = s1.configuration.box[0]
    ly1 = s1.configuration.box[1]
    lx2 = s2.configuration.box[0]
    ly2 = s2.configuration.box[1]
    if lx1 != lx2:
        print("Error, Lx=%g,%g" % (lx1, lx2))
        sys.exit(1)
    n1 = s1.particles.N
    N = n1 + s2.particles.N
    Lx, Ly = lx1, ly1 + ly2
    pos = np.zeros((N, 3), dtype=np.float32)
    type_id = np.zeros(N, dtype=np.uint32)

    dy1 = -ly2/2
    dy2 = ly1/2
    pos[:n1, 0] = s1.particles.position[:, 0]
    pos[n1:, 0] = s2.particles.position[:, 0]
    pos[:n1, 1] = s1.particles.position[:, 1] + dy1
    pos[n1:, 1] = s2.particles.position[:, 1] + dy2
    pos[:n1, 2] = s1.particles.position[:, 2]
    pos[n1:, 2] = s2.particles.position[:, 2]
    type_id[:n1] = s1.particles.typeid
    type_id[n1:] = s2.particles.typeid


    s2 = hoomd.Snapshot()
    s2.configuration.box = [Lx, Ly, 1, 0, 0, 0]
    s2.particles.N = N
    s2.particles.position = pos
    s2.particles.typeid = type_id
    s2.particles.types = s1.particles.types
    s2.configuration.step = 0
    return s2

def rot90(s: hoomd.Snapshot) -> hoomd.Snapshot:
    # Ly = s.configuration.box[1]
    # Lx = s.configuration.box[0]

    x = s.particles.position[:, 0].copy()
    y = s.particles.position[:, 1].copy()
    theta = s.particles.position[:, 2]

    # x, y = -y, x
    theta += np.pi / 2
    theta[theta>np.pi] -= np.pi * 2
    s.particles.position[:, 0] = -y
    s.particles.position[:, 1] = x
    s.particles.position[:, 2] = theta
    return s


def split(s: hoomd.Snapshot, x_shift=0):
    s = shift_pos(s, x_shift, 0)

    Lx = s.configuration.box[0]
    Ly = s.configuration.box[1]

    x = s.particles.position[:, 0]
    mask_left = x < 0
    mask_right = x >= 0

    n_left = np.sum(mask_left)
    n_right = np.sum(mask_right)

    print(n_left, n_right, n_left+n_right - x.size)

    new_Lx = Lx / 2
    new_area = new_Lx * Ly
    
    s_left = hoomd.Snapshot()
    s_left.configuration.box = [new_Lx, Ly, 1, 0, 0, 0]
    s_left.particles.N = n_left
    s_left.particles.position = s.particles.position[mask_left]
    s_left.particles.typeid = s.particles.typeid[mask_left]
    s_left.particles.types = s.particles.types
    s_left.configuration.step = 0
    s_left.particles.position[:, 0] += new_Lx / 2

    s_right = hoomd.Snapshot()
    s_right.configuration.box = [new_Lx, Ly, 1, 0, 0, 0]
    s_right.particles.N = n_right
    s_right.particles.position = s.particles.position[mask_right]
    s_right.particles.typeid = s.particles.typeid[mask_right]
    s_right.particles.types = s.particles.types
    s_right.configuration.step = 0
    s_right.particles.position[:, 0] -= new_Lx / 2

    nA = np.sum(s.particles.typeid == 0)
    nB = np.sum(s.particles.typeid == 1)

    nA_left = np.sum(s_left.particles.typeid == 0)
    nB_left = np.sum(s_left.particles.typeid == 1)
    print(nA_left + nB_left, n_left)

    rhoA_left = np.round(nA_left / new_area * 8, 2) / 8
    rhoB_left = np.round(nB_left / new_area * 8, 2) / 8

    nA_left_new = np.round(rhoA_left * new_area)
    nB_left_new = np.round(rhoB_left * new_area)

    print(rhoA_left, nA_left_new, nA_left, rhoB_left, nB_left_new, nB_left)

    nA_left_new = int(nA_left_new)
    nB_left_new = int(nB_left_new)

    print((nA_left_new - nA_left)/nA_left)
    print((nB_left_new - nB_left)/nB_left)


    nA_right = np.sum(s_right.particles.typeid == 0)
    nB_right = np.sum(s_right.particles.typeid == 1)
    rhoA_right = np.round(nA_right / new_area * 8, 2) / 8
    rhoB_right = np.round(nB_right / new_area * 8, 2) / 8

    nA_right_new = rhoA_right * new_area
    nB_right_new = rhoB_right * new_area

    print(rhoA_right, nA_right_new, rhoB_right, nB_right_new)

    nA_right_new = int(nA_right_new)
    nB_right_new = int(nB_right_new)
    print((nA_right_new - nA_right)/nA_right)
    print((nB_right_new - nB_right)/nB_right)

    print(nA_left_new - nA_left)
    print(nB_left_new - nB_left)
    print(nA_right_new - nA_right)
    print(nB_right_new - nB_right)

    s_left = adjust_density(s_left, rhoA_left, rhoB_left, mode="copy")
    s_right = adjust_density(s_right, rhoA_right, rhoB_right, mode="copy")

    rhoA_left = np.sum(s_left.particles.typeid == 0) / new_area
    rhoB_left = np.sum(s_left.particles.typeid == 1) / new_area
    rhoA_right = np.sum(s_right.particles.typeid == 0) / new_area
    rhoB_right = np.sum(s_right.particles.typeid == 1) / new_area

    print(rhoA_left, rhoB_left)
    print(rhoA_right, rhoB_right)
    return s_left, s_right


if __name__ == "__main__":
    folder = "/scratch03.local/yduan/topoVM/dissenters/L800_new2/"
    # folder = "build/data"
    basename = "L800_800_d0.1000_e0.300_r1_s2011.gsd"

    fname = f"{folder}/{basename}"
    snap = read_one_frame(fname, -1)
    # snap = duplicate(snap, 2, 2)
    # snap = scale(snap, 2, 2, 0.25)
    snap = add_A_particles(snap, nA_added=32000)
    snap.configuration.step = 0

    snap = sort_by_type(snap)
    fout = f"{folder}/L800_800_d0.095238_e0.300_r1.05_s2001.gsd"
    f = hoomd.open(name=fout, mode='wb')
    f.append(snap)

