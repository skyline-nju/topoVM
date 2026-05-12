#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "rand.h"
#include "comn.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Vector_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Vector_2<K> Vec_2;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;
typedef CGAL::Periodic_2_triangulation_vertex_base_2<GT> Vb;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, GT, Vb> VbInfo;
typedef CGAL::Periodic_2_triangulation_face_base_2<GT> Fb;
typedef CGAL::Triangulation_data_structure_2<VbInfo, Fb> Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT, Tds> PDT;
typedef PDT::Point Point;
typedef Tds::Vertex_handle Vertex_handle;

/****************************************************************************
 * Binary mixture of self-propelled particles aligning with Voronoi neighbors
 * The position and orientation of particle i from species S are updated as:
 * x_i(t+h) = x_i(t) + v0 * h * u_i(t),
 * theta_i(t+h) = sum_{j in Voronoi neighbors} J_{SS'} sin[theta_j(t) - theta_i(t)] * h + sqrt(24 * Dr * h) * xi,
 * where S' is the species of particle j, 
 * J_{SS'} is the strength of aligning interaction felt by S particles from S' particles, 
 * xi is a random number uniformly distributed in [-0.5, 0.5],
 * Dr is the rotational diffusion coefficient,
 * h is the time step.
 * 
 ****************************************************************************/
class VM_Voro_Conti_2 {
public:
  VM_Voro_Conti_2(double Lx, double Ly,
      double J_AA, double J_AB, double J_BA, double J_BB,
      double Dr, double h, double phi_A, double v0 = 1);

  ~VM_Voro_Conti_2() { delete DT_; }

  int get_n() const { return n_; }
  int get_nA() const { return nA_; }

  template <typename TPos, typename TVel>
  void ini_pos_vel(const TPos *x, const TPos *y, const TVel *vx, const TVel *vy);

  template <typename T>
  void get_x(int i, T &x, T &y) const { auto p = vh_arr_[i]->point(); x = p.x(); y = p.y(); }

  template <typename T>
  void get_theta(int i, T &theta) const { theta = atan2(v_arr_[i].y, v_arr_[i].x);}

  unsigned int get_type(int i) const { return vh_arr_[i]->info() < nA_ ? 0 : 1; }

  void cal_order_para(double &phi_A, double &phi_B, double &theta_A, double &theta_B) const;

  template <typename T>
  void get_pos_arr(T *pos, uint32_t *type_id = nullptr) const;

  template <typename TRan>
  void ini_rand(TRan &myran);

  template <typename TRan>
  void ini_rand(TRan &myran, double theta0_A, double theta0_B);

  template <typename TSnap>
  void ini_from_snap(TSnap& snap_reader);

  template <typename TSnap>
  void dump(int i_step, TSnap& snap_writer);

  void pairwise_align(unsigned int i, unsigned int j);

  void collide();

  template <typename TRan>
  void update_vel_pos(unsigned int i, TRan &myran, double &x, double &y);


  template<typename TRan>
  void update_vertex_additively(TRan &myran);

  template<typename TRan>
  void update_vertex(TRan &myran);

  template<typename TRan>
  void update_vertex_w_spatial_sort(TRan &myran);

  template<typename TRan>
  void update_vertex_by_moving(TRan &myran);

  template <typename TRan>
  void stream(TRan &myran) {
    update_vertex_additively(myran); 
    // update_vertex(myran);
    // update_vertex_w_spatial_sort(myran);
    // update_vertex_by_moving(myran);
  }


protected:
  double v0_h_;
  double Jh_[4]{};      // J_SS' * h, stored in the order of AA, AB, BA and BB
  double sqrt_24_Dr_h_; // sqrt(24 * Dr * h)
  double Lx_;
  double Ly_;
  int nA_;
  int n_;

  // array of velocity
  std::vector<float2> v_arr_;
  std::vector<float> tau_arr_;

  // array of vertex handle
  std::vector<Vertex_handle> vh_arr_;

  PDT *DT_;
};

template <typename TPos, typename TVel>
void VM_Voro_Conti_2::ini_pos_vel(const TPos *x, const TPos *y, const TVel *vx, const TVel *vy) {
  v_arr_.reserve(n_);
  tau_arr_.reserve(n_);
  vh_arr_.reserve(n_);

  for (int i = 0; i < n_; i++) {
      v_arr_.emplace_back(vx[i], vy[i]);
      tau_arr_.push_back(0);
      auto vh = DT_->insert(Point(x[i], y[i]));
      vh->info() = i;
      vh_arr_.push_back(vh);
  }

}


template <typename T>
void VM_Voro_Conti_2::get_pos_arr(T *pos, uint32_t *type_id) const {
  for (int i = 0; i < n_; i++) {
    size_t j = i * 3;
    double x, y, theta;
    get_x(i, x, y);
    get_theta(i, theta);
    pos[j] = x - Lx_ / 2;
    pos[j+1] = y - Ly_ / 2;
    pos[j+2] = theta;
    if (type_id) {
        type_id[i] = get_type(i);
    }
  }
}

template <typename TRan>
void VM_Voro_Conti_2::ini_rand(TRan &myran) {
  double *x = new double[n_];
  double *y = new double[n_];
  double *vx = new double[n_];
  double *vy = new double[n_];
  for (int i = 0; i < n_; i++) {
      x[i] = myran.doub() * Lx_;
      y[i] = myran.doub() * Ly_;
      double theta = myran.doub() * 2 * PI;
      vx[i] = cos(theta);
      vy[i] = sin(theta);
  }
  ini_pos_vel(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <typename TRan>
void VM_Voro_Conti_2::ini_rand(TRan &myran, double theta0_A, double theta0_B) {
  double *x = new double[n_];
  double *y = new double[n_];
  double *vx = new double[n_];
  double *vy = new double[n_];
  double vx_A = cos(theta0_A);
  double vy_A = sin(theta0_A);
  double vx_B = cos(theta0_B);
  double vy_B = sin(theta0_B);
  for (int i = 0; i < n_; i++) {
    x[i] = myran.doub() * Lx_;
    y[i] = myran.doub() * Ly_;
    if (i < nA_) {
        vx[i] = vx_A;
        vy[i] = vy_A;
    } else {
        vx[i] = vx_B;
        vy[i] = vy_B;
    }
  }
  ini_pos_vel(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <typename TSnap>
void VM_Voro_Conti_2::ini_from_snap(TSnap &snap_reader) {
  double *x = new double[n_];
  double *y = new double[n_];
  double *vx = new double[n_];
  double *vy = new double[n_];
  snap_reader.read_last_frame(x, y, vx, vy);
  ini_pos_vel(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <typename TSnap>
void VM_Voro_Conti_2::dump(int i_step, TSnap &snap_writer) {
  if (snap_writer.need_export(i_step)) {
    float* pos = new float[n_*3];
    uint32_t *type_id = new uint32_t[n_];
    get_pos_arr(pos, type_id);
    snap_writer.dump(i_step, n_, pos, type_id);
    delete [] pos;
    delete [] type_id;
  }
}


template <typename TRan>
void VM_Voro_Conti_2::update_vel_pos(unsigned int i, TRan &myran, double &x, double &y) {
  // update velocity of particle i
  double noise = (myran.doub() - 0.5) * sqrt_24_Dr_h_;
  double d_theta = tau_arr_[i] + noise;
  double s = sin(d_theta);
  double c = cos(d_theta);
  double vx = v_arr_[i].x * c - v_arr_[i].y * s;
  double vy = v_arr_[i].x * s + v_arr_[i].y * c;
  double tmp = 1.0 / sqrt(vx * vx + vy * vy);
  v_arr_[i].x = vx * tmp;
  v_arr_[i].y = vy * tmp;
  tau_arr_[i] = 0;

  // update position of particle i
  get_x(i, x, y);
  x += v_arr_[i].x * v0_h_;
  y += v_arr_[i].y * v0_h_;
  if (x < 0) {
      x += Lx_;
  } else if (x >= Lx_) {
      x -= Lx_;
  }
  if (y < 0) {
      y += Ly_;
  } else if (y >= Ly_) {
      y -= Ly_;
  }

}

template <typename TRan>
void VM_Voro_Conti_2::update_vertex_additively(TRan &myran) {
  for (int i = 0; i < n_; i++){
      double x, y;
      update_vel_pos(i, myran, x, y);
      Point p_new(x, y);
      auto vh_old = vh_arr_[i];
      vh_arr_[i] = DT_->insert(p_new, vh_old->face());
      vh_arr_[i]->info() = vh_old->info();
      DT_->remove(vh_old);
  }
}

template <typename TRan>
void VM_Voro_Conti_2::update_vertex(TRan &myran) {
  // 1. 先更新所有粒子的位置和速度（不修改三角剖分）
  std::vector<Point> new_pos(n_);
  for (int i = 0; i < n_; i++) {
    double x, y;
    update_vel_pos(i, myran, x, y);
    new_pos[i] = Point(x, y);
  }
  // 2. 清空旧的三角剖分，批量插入新位置
  DT_->clear();
  vh_arr_.clear();
  for (int i = 0; i < n_; i++) {
    auto vh = DT_->insert(new_pos[i]);
    vh->info() = i;
    vh_arr_.push_back(vh);
  }
}

template <typename TRan>
void VM_Voro_Conti_2::update_vertex_w_spatial_sort(TRan &myran) {
   // 1. 先更新所有粒子的位置和速度（不修改三角剖分）
  std::vector<std::pair<Point, unsigned int>> point_info_pairs;
  point_info_pairs.reserve(n_);
  for (int i = 0; i < n_; i++) {
    double x, y;
    update_vel_pos(i, myran, x, y);
    Point p(x, y);
    point_info_pairs.push_back(std::make_pair(p, i));
  }

  // 2. 清空旧的三角剖分，批量插入新位置
  DT_->clear();
  DT_->insert(point_info_pairs.begin(), point_info_pairs.end(), true);
  for (auto vit = DT_->finite_vertices_begin(); vit != DT_->finite_vertices_end(); ++vit) {
       unsigned int original_idx = vit->info();
       vh_arr_[original_idx] = vit;
  }

  // for (int i = 0; i < n_; i++) {
  //   assert(vh_arr_[i]->info() == i); // 确保原始索引匹配
  // }
}

template <typename TRan>
void VM_Voro_Conti_2::update_vertex_by_moving(TRan &myran) {
  for (int i = 0; i < n_; i++){
      double x, y;
      update_vel_pos(i, myran, x, y);
      Point p_new(x, y);
      vh_arr_[i] = DT_->move_point(vh_arr_[i], p_new);
      // auto vh_old = vh_arr_[i];
      // vh_arr_[i] = DT_->insert(p_new);
      vh_arr_[i]->info() = i;
      // DT_->remove(vh_old);
      if (vh_arr_[i]->info() != i) {
        std::cout << "Error in update_vertex_by_moving, vertex info does not match with i!" << std::endl;
        std::cout << vh_arr_[i]->info() << " != " << i << std::endl;
        exit(1);
      }
  }
}
