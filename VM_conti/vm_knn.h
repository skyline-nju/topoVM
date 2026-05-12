#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include "rand.h"
#include "comn.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;
typedef boost::tuple<Point_2, int> Point_and_int;
typedef CGAL::Search_traits_2<Kernel> Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int,
                                    CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
                                    Traits_base> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;
typedef K_neighbor_search::Distance Distance;


class VM_kNN {
public:
  VM_kNN(double Lx, double Ly, double rho0, double J0, double Dr, double h, double v0=1.0);

  int get_n() const { return N_; }

  template <typename T>
  void get_x(int i, T &x, T &y) const { x = pos_arr_[i].x(); y = pos_arr_[i].y(); }

  template <typename T>
  void get_theta(int i, T &theta) const { theta = atan2(v_arr_[i].y, v_arr_[i].x);}

  template <typename TPos, typename TVel>
  void ini_pos_vel(const TPos *x, const TPos *y, const TVel *vx, const TVel *vy);

  template <typename T>
  void get_pos_arr(T *pos) const;

  template <typename TRan>
  void ini_rand(TRan &myran);

  template <typename TRan>
  void ini_rand(TRan &myran, double theta0);

  template <typename TSnap>
  void ini_from_snap(TSnap& snap_reader);

  template <typename TSnap>
  void dump(int i_step, TSnap& snap_writer);

  void add_padded_particles();

  void add_padded_particles(double dx, double dy);

  void del_padded_particles() {pos_arr_.resize(N_); idx_arr_.resize(N_);}

  void align_one_side(unsigned int i, unsigned int j) {
      // J0 * sin (theta_j - theta_i) * h
      tau_arr_[i] += J0_h_ * (v_arr_[j].y * v_arr_[i].x - v_arr_[j].x * v_arr_[i].y); 
  }

  void align_base(int k, double *pt_max_dis2=nullptr);

  void align(int k) { add_padded_particles(); align_base(k); del_padded_particles(); }

  void align_padded(int k);

  void collide(int k) { align_padded(k); }

  template <typename TRan>
  void stream(TRan& myran);

  void cal_order_para(double &phi, double &theta) const;

protected:
  double Lx_;
  double Ly_;
  int N_;
  double J0_h_;
  double v0_h_;
  double sqrt_24_Dr_h_; // sqrt(24 * Dr * h)


  // array of velocity
  std::vector<float2> v_arr_;
  std::vector<float> tau_arr_;

  // array of pos
  std::vector<Point_2> pos_arr_;
  std::vector<int> idx_arr_; // the index of particle corresponding to each point in pos_arr_

  double padded_dx2_  = 0;
};

template <typename TPos, typename TVel>
void VM_kNN::ini_pos_vel(const TPos *x, const TPos *y, const TVel *vx, const TVel *vy) {
  for (int i = 0; i < N_; i++) {
    v_arr_.emplace_back(vx[i], vy[i]);
    pos_arr_.emplace_back(x[i], y[i]);
    idx_arr_.push_back(i);
  }
}


template <typename T>
void VM_kNN::get_pos_arr(T *pos) const {
   for (int i = 0; i < N_; i++) {
    size_t j = i * 3;
    double x, y, theta;
    get_x(i, x, y);
    get_theta(i, theta);
    pos[j  ] = x - 0.5 * Lx_;
    pos[j+1] = y - 0.5 * Ly_;
    pos[j+2] = theta;
  }
}

template <typename TRan>
void VM_kNN::ini_rand(TRan &myran) {
  double *x = new double[N_];
  double *y = new double[N_];
  double *vx = new double[N_];
  double *vy = new double[N_];
  for (int i = 0; i < N_; i++) {
    x[i] = myran.doub() * Lx_;
    y[i] = myran.doub() * Ly_;
    double theta = PI * myran.doub() * 2.;
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
void VM_kNN::ini_rand(TRan &myran, double theta0) {
  double *x = new double[N_];
  double *y = new double[N_];
  double *vx = new double[N_];
  double *vy = new double[N_];
  double vx0 = cos(theta0);
  double vy0 = sin(theta0);
  for (int i = 0; i < N_; i++) {
    x[i] = myran.doub() * Lx_;
    y[i] = myran.doub() * Ly_;
    vx[i] = vx0;
    vy[i] = vy0;
  }
  ini_pos_vel(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <typename TSnap>
void VM_kNN::ini_from_snap(TSnap &snap_reader) {
  double *x = new double[N_];
  double *y = new double[N_];
  double *vx = new double[N_];
  double *vy = new double[N_];
  snap_reader.read_last_frame(x, y, vx, vy);
  ini_pos_vel(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <typename TSnap>
void VM_kNN::dump(int i_step, TSnap &snap_writer) {
  if (snap_writer.need_export(i_step)) {
    float* pos = new float[N_*3];
    get_pos_arr(pos);
    snap_writer.dump(i_step, N_, pos);
    delete [] pos;
  }
}

template <typename TRan>
void VM_kNN::stream(TRan &myran) {
  for (int i = 0; i < N_; i++) {
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
    double x, y;
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
    pos_arr_[i] = Point_2(x, y);
  }
}