#include "vm_knn.h"

VM_kNN::VM_kNN(double Lx, double Ly, double rho0, double J0, double Dr, double h, double v0)
  : Lx_(Lx), Ly_(Ly), J0_h_(J0 * h), v0_h_(v0 * h), sqrt_24_Dr_h_(sqrt(24.0 * Dr * h)) {
    N_ = static_cast<int>(rho0 * Lx_ * Ly_);
    v_arr_.reserve(N_);
    tau_arr_.reserve(N_);
    pos_arr_.reserve(N_ * 9);
    idx_arr_.reserve(N_ * 9);
}


void VM_kNN::add_padded_particles() {
  double dx_arr[8] = {-Lx_, -Lx_, -Lx_, 0, 0, Lx_, Lx_, Lx_};
  double dy_arr[8] = {-Ly_, 0, Ly_, -Ly_, Ly_, -Ly_, 0, Ly_};

  for (int i = 0; i < N_; i++) {
    double x, y;
    get_x(i, x, y);
    for (int j = 0; j < 8; j++) {
      pos_arr_.emplace_back(x + dx_arr[j], y + dy_arr[j]);
      idx_arr_.emplace_back(i);
    }
  }
}

void VM_kNN::add_padded_particles(double dx, double dy) {
  for (int i = 0; i < N_; i++) {
    double x, y;
    get_x(i, x, y);

    if (x <= dx) {
      pos_arr_.emplace_back(x + Lx_, y);
      idx_arr_.emplace_back(i);
      if (y <= dy) {
        pos_arr_.emplace_back(x, y+Ly_);
        idx_arr_.emplace_back(i);
        pos_arr_.emplace_back(x + Lx_, y + Ly_);
        idx_arr_.emplace_back(i);
      } else if (y > Ly_ - dy) {
        pos_arr_.emplace_back(x, y-Ly_);
        idx_arr_.emplace_back(i);
        pos_arr_.emplace_back(x + Lx_, y - Ly_);
        idx_arr_.emplace_back(i);
      }
    } else if (x <= Lx_ - dx) {
      if (y <= dy) {
        pos_arr_.emplace_back(x, y + Ly_);
        idx_arr_.emplace_back(i);
      } else if (y > Ly_ - dy) {
        pos_arr_.emplace_back(x, y - Ly_);
        idx_arr_.emplace_back(i);
      }
    } else {
      pos_arr_.emplace_back(x - Lx_, y);
      idx_arr_.emplace_back(i);
      if (y <= dy) {
        pos_arr_.emplace_back(x, y + Ly_);
        idx_arr_.emplace_back(i);
        pos_arr_.emplace_back(x - Lx_, y + Ly_);
        idx_arr_.emplace_back(i);
      } else if (y > Ly_ - dy) {
        pos_arr_.emplace_back(x, y - Ly_);
        idx_arr_.emplace_back(i);
        pos_arr_.emplace_back(x - Lx_, y - Ly_);
        idx_arr_.emplace_back(i);
      }
    }
  }
}


void VM_kNN::align_base(int k, double *pt_max_dis2) {
  Tree tree(boost::make_zip_iterator(boost::make_tuple( pos_arr_.begin(),idx_arr_.begin())),
          boost::make_zip_iterator(boost::make_tuple( pos_arr_.end(),idx_arr_.end())));

  for (int i = 0; i < N_; i++) {
    K_neighbor_search search(tree, pos_arr_[i], k+1);
    for(K_neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
      if (it != search.begin()) {  // The first neighbor is particle i itself and should be excluded 
        int j = boost::get<1>(it->first);
        align_one_side(i, j);
        if (pt_max_dis2 && *pt_max_dis2 < it->second) {
          *pt_max_dis2 = it->second;
        }
      }
    }
  }
}

void VM_kNN::align_padded(int k) {
  if (padded_dx2_ == 0 || padded_dx2_ >= Lx_ * Lx_ || padded_dx2_ >= Ly_ * Ly_) {
    add_padded_particles();
    padded_dx2_ = 0;
    align_base(k, &padded_dx2_);
    del_padded_particles();
  } else {
    double dx = sqrt(padded_dx2_) + 2 * v0_h_;
    double dy = dx;
    add_padded_particles(dx, dy);
    padded_dx2_ = 0;
    align_base(k, &padded_dx2_);
    del_padded_particles();
  }
}

void VM_kNN::cal_order_para(double &phi, double &theta) const {
  double vx = 0.;
  double vy = 0.;
  for (int i = 0; i < N_; i++) {
    vx += v_arr_[i].x;
    vy += v_arr_[i].y;
  }
  phi = sqrt(vx * vx + vy * vy) / N_;
  theta = atan2(vy, vx);
}
