#include "vm_voro_2.h"
#include <CGAL/spatial_sort.h>

struct Extract_point {
  typedef Point result_type;
  const Point& operator()(const std::pair<Point, unsigned int>& p) const {
    return p.first;
  }
};

VM_Voro_Conti_2::VM_Voro_Conti_2(double Lx, double Ly,
                                 double J_AA, double J_AB, double J_BA, double J_BB,
                                 double Dr, double h,
                                 double phi_A, double v0) 
  : Lx_(Lx), Ly_(Ly), v0_h_(h * v0) {
  Jh_[0] = J_AA * h;
  Jh_[1] = J_AB * h;
  Jh_[2] = J_BA * h;
  Jh_[3] = J_BB * h;
  sqrt_24_Dr_h_ = std::sqrt(24 * Dr * h);
  n_ = static_cast<int>(Lx * Ly);
  nA_ = static_cast<int>(phi_A * Lx * Ly);

  PDT::Iso_rectangle domain(0, 0, Lx_, Ly_);
  DT_ = new PDT(domain);
}

void VM_Voro_Conti_2::cal_order_para(double &phi_A, double &phi_B,
                                     double &theta_A, double &theta_B) const {
  double vx_A = 0.;
  double vy_A = 0.;
  double vx_B = 0.;
  double vy_B = 0.;
  for (int i = 0; i < nA_; i++) {
    vx_A += v_arr_[i].x;
    vy_A += v_arr_[i].y;
  }
  phi_A = sqrt(vx_A * vx_A + vy_A * vy_A) / nA_;
  theta_A = atan2(vy_A, vx_A);

  for (int i = nA_; i < n_; i++) {
    vx_B += v_arr_[i].x;
    vy_B += v_arr_[i].y;
  }
  if (n_ - nA_ > 0) {
    phi_B = sqrt(vx_B * vx_B + vy_B * vy_B) / (n_ - nA_);
    theta_B = atan2(vy_B, vx_B); 
  } else {
    phi_B = 0;
    theta_B = 0;
  }


}

void VM_Voro_Conti_2::pairwise_align(unsigned int i, unsigned int j) {
  double sin_dtheta = v_arr_[j].y * v_arr_[i].x - v_arr_[j].x * v_arr_[i].y; // sin (theta_j - theta_i)
  int s1 = (i < nA_) ? 0 : 1;
  int s2 = (j < nA_) ? 0 : 1;
  // if (s1 != get_type(i) || s2 != get_type(j)) {
  //   std::cout << "Error in pairwise_align, type_id does not match with i and j!" << std::endl;
  //   exit(1);
  // }
  double J12_h = Jh_[s1 * 2 + s2];
  double J21_h = Jh_[s2 * 2 + s1];
  tau_arr_[i] += J12_h * sin_dtheta;
  tau_arr_[j] -= J21_h * sin_dtheta;
}


void VM_Voro_Conti_2::collide() {
  for (auto fit = DT_->periodic_triangles_begin(PDT::UNIQUE);
            fit != DT_->periodic_triangles_end(PDT::UNIQUE); ++fit)
  {
    unsigned int idx0 = fit.get_face()->vertex(0)->info();
    unsigned int idx1 = fit.get_face()->vertex(1)->info();
    unsigned int idx2 = fit.get_face()->vertex(2)->info();
    if (idx0 < idx1) {
        pairwise_align(idx0, idx1);
    }
    if (idx1 < idx2) {
        pairwise_align(idx1, idx2);
    }
    if (idx2 < idx0) {
        pairwise_align(idx2, idx0);
    }
  }
}


