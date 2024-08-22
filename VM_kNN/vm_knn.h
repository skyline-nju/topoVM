#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <utility>

#include "particle.h"
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

// Base class for Vicsek model
class VM {
public:
  VM(double Lx, double Ly, int N, double eta, double v0=0.5)
  : Lx_(Lx), Ly_(Ly), half_Lx_(Lx * 0.5), half_Ly_(Ly * 0.5), N_(N), eta_(eta), v0_(v0) {}

  ~VM(){};
 
  // get data
  int get_n_birds() const { return N_; }

protected:
  double v0_;
  double Lx_;
  double Ly_;
  double half_Lx_;
  double half_Ly_;
  int N_;
  double eta_;
};

/*************************************************************************************
 *      Class for one-species flocks with kNN interaction
 * 
 * 
 ************************************************************************************/
template <class BaseV>
class VM_kNN: public VM {
public:
  VM_kNN(double Lx, double Ly, int N, double eta, double v0=0.5);

  // ~VM_kNN(){}

  template <typename T>
  void input_data(const T *x, const T *y, const T *vx, const T *vy);

  template <typename T>
  void get_pos_arr(T* pos) const;

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

  void align_base(int k, double *pt_max_dis2=nullptr);

  void align(int k);

  void align_padded(int k);

  template<typename TRan>
  void stream(TRan& myran);

  void get_x(int i, double &x, double &y) const;
  void get_v(int i, double &vx, double &vy) const { v_arr_[i].get_v(vx, vy); }
  void get_theta(int i, double & Theta) const { v_arr_[i].get_theta(Theta); }
  template <typename T>
  void get_order_para(T &phi, T &theta) const;

protected:
  // array of velocity
  std::vector<BaseV> v_arr_;

  // array of pos
  std::vector<Point_2> pos_arr_;

  // array of idx
  std::vector<int> idx_arr_;

  double padded_dx2_ = 0;
};

template <class BaseV>
VM_kNN<BaseV>::VM_kNN(double Lx, double Ly, int N, double eta, double v0)
  : VM(Lx, Ly, N, eta, v0) {
    v_arr_.reserve(N_);

    size_t N_new = N_ * 9;
    pos_arr_.reserve(N_new);
    idx_arr_.reserve(N_new);
}

template<class BaseV>
template<typename T>
void VM_kNN<BaseV>::input_data(const T *x, const T *y, const T *vx, const T *vy) {
  for (int i = 0; i < N_; i++) {
    v_arr_.emplace_back(vx[i], vy[i]);
    pos_arr_.emplace_back(x[i], y[i]);
    idx_arr_.push_back(i);
  }
}

template <class BaseV>
template <typename T>
void VM_kNN<BaseV>::get_pos_arr(T *pos) const {
  for (int i = 0; i < N_; i++) {
    size_t j = i * 3;
    double x, y, theta;
    get_x(i, x, y);
    get_theta(i, theta);
    pos[j] = x - half_Lx_;
    pos[j+1] = y - half_Ly_;
    pos[j+2] = theta;
  }
}

template <class BaseV>
template <typename TRan>
void VM_kNN<BaseV>::ini_rand(TRan &myran) {
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
  input_data(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <class BaseV>
template <typename TRan>
void VM_kNN<BaseV>::ini_rand(TRan &myran, double theta0) {
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
  input_data(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <class BaseV>
template <typename TSnap>
void VM_kNN<BaseV>::ini_from_snap(TSnap &snap_reader) {
  double *x = new double[N_];
  double *y = new double[N_];
  double *vx = new double[N_];
  double *vy = new double[N_];
  snap_reader.read_last_frame(x, y, vx, vy);
  input_data(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <class BaseV>
template <typename TSnap>
void VM_kNN<BaseV>::dump(int i_step, TSnap &snap_writer) {
  if (snap_writer.need_export(i_step)) {
    float* pos = new float[N_*3];
    get_pos_arr(pos);
    snap_writer.dump(i_step, N_, pos);
    delete [] pos;
  }
}

template <class BaseV>
void VM_kNN<BaseV>::add_padded_particles() {
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


template <class BaseV>
void VM_kNN<BaseV>::add_padded_particles(double dx, double dy) {
  for (int i = 0; i < N_; i++) {
    double x, y;
    get_x(i, x, y);

    if (x <= dx) {
      pos_arr_.emplace_back(x + Lx_, y);
      idx_arr_.emplace_back(i);
      if (y <= dy) {
        pos_arr_.emplace_back(x + Lx_, y + Ly_);
        idx_arr_.emplace_back(i);
      } else if (y > Ly_ - dy) {
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
        pos_arr_.emplace_back(x - Lx_, y + Ly_);
        idx_arr_.emplace_back(i);
      } else if (y > Ly_ - dy) {
        pos_arr_.emplace_back(x - Lx_, y - Ly_);
        idx_arr_.emplace_back(i);
      }
    }
  }

}

template <class BaseV>
void VM_kNN<BaseV>::align_base(int k, double *pt_max_dis2) {
    Tree tree(boost::make_zip_iterator(boost::make_tuple( pos_arr_.begin(),idx_arr_.begin())),
            boost::make_zip_iterator(boost::make_tuple( pos_arr_.end(),idx_arr_.end())));
  
    for (int i = 0; i < N_; i++) {
      K_neighbor_search search(tree, pos_arr_[i], k);
      for(K_neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
        int j = boost::get<1>(it->first);
        v_arr_[i].collide_one_side(v_arr_[j]);
        if (pt_max_dis2 && *pt_max_dis2 < it->second) {
          *pt_max_dis2 = it->second;
        }
      }
    }
}

template <class BaseV>
void VM_kNN<BaseV>::align(int k) {
  add_padded_particles();
  align_base(k);
  del_padded_particles();
}


template <class BaseV>
void VM_kNN<BaseV>::align_padded(int k) {
  if (padded_dx2_ == 0) {
    add_padded_particles();
    align_base(k, &padded_dx2_);
    del_padded_particles();
  } else {
    double dx = sqrt(padded_dx2_) + 1;
    double dy = dx;
    add_padded_particles(dx, dy);
    padded_dx2_ = 0;
    align_base(k, &padded_dx2_);
    del_padded_particles();
  }
}
template <class BaseV>
template <typename TRan>
void VM_kNN<BaseV>::stream(TRan &myran) {
  for (int i = 0; i < N_; i++) {
    v_arr_[i].update_v(eta_, myran);
    double x, y;
    get_x(i, x, y);
    v_arr_[i].update_x(x, y, v0_, Lx_, Ly_);
    pos_arr_[i] = Point_2(x, y);
  }
}

template<class BaseV>
void VM_kNN<BaseV>::get_x(int i, double & x, double & y) const {
  x = pos_arr_[i].x();
  y = pos_arr_[i].y();
}

template <class BaseV>
template <typename T>
void VM_kNN<BaseV>::get_order_para(T &phi, T &theta) const {
  double vx = 0;
  double vy = 0;
  for (int i = 0; i < N_; i++) {
    vx += v_arr_[i].vx;
    vy += v_arr_[i].vy;
  }
  vx /= N_;
  vy /= N_;
  phi = sqrt(vx * vx + vy * vy);
  theta = atan2(vy, vx); 
}

/***************************************************************************************
 *      Class for binary mixture of aligners and dissenters with voronoi-type 
 *      interaction
 * 
 ***************************************************************************************/
/*
template <class BaseV>
class VM_Voro_AlignerDissenter: public VM_Voro<BaseV> {
public:
  VM_Voro_AlignerDissenter(double Lx, double Ly, int N, double eta, double v0, int n_dis)
    : VM_Voro<BaseV>(Lx, Ly, N, eta, v0), n_dis_(n_dis) {}


  template <typename TRan>
  void ini_rand(TRan &myran) {VM_Voro<BaseV>::ini_rand(myran);}
  template <typename TRan>
  void ini_rand(TRan &myran, double theta0);
  template <typename TSnap>
  void ini_from_snap(TSnap& snap_reader);
  
  void align();

  void align_nearest_neighbor();

  template <typename T>
  void get_order_para(T &phi, T &theta) const;

  template <typename TSnap>
  void dump(int i_step, TSnap& snap_writer);
protected:
  int n_dis_;  // the first n_dis_ particles out of N_ particles are dissenters, 
               // while the rest is aligners
};

template <class BaseV>
void VM_Voro_AlignerDissenter<BaseV>::align() {
  for (auto fit = VM_Voro<BaseV>::DT_->periodic_triangles_begin(PDT::UNIQUE);
    fit != VM_Voro<BaseV>::DT_->periodic_triangles_end(PDT::UNIQUE); ++fit) {
    unsigned int idx0 = fit.get_face()->vertex(0)->info();
    unsigned int idx1 = fit.get_face()->vertex(1)->info();
    unsigned int idx2 = fit.get_face()->vertex(2)->info();
    bool is_aligner_0 = idx0 >= n_dis_;
    bool is_aligner_1 = idx1 >= n_dis_;
    bool is_aligner_2 = idx2 >= n_dis_;
    if (idx0 < idx1) {
      VM_Voro<BaseV>::v_arr_[idx0].collide(&VM_Voro<BaseV>::v_arr_[idx1], is_aligner_0, is_aligner_1);
    }
    if (idx1 < idx2) {
      VM_Voro<BaseV>::v_arr_[idx1].collide(&VM_Voro<BaseV>::v_arr_[idx2], is_aligner_1, is_aligner_2);
    }
    if (idx2 < idx0) {
      VM_Voro<BaseV>::v_arr_[idx2].collide(&VM_Voro<BaseV>::v_arr_[idx0], is_aligner_2, is_aligner_0);
    }
  }
}

// TODO
template <class BaseV>
void VM_Voro_AlignerDissenter<BaseV>::align_nearest_neighbor() {
  for (auto pit = VM_Voro<BaseV>::DT_->periodic_points_begin(PDT::UNIQUE);
    pit != VM_Voro<BaseV>::DT_->periodic_points_end(PDT::UNIQUE); ++pit) {
    unsigned int idx0 = pit.get_vertex()->info();
    bool is_aligner_0 = idx0 >= n_dis_;
    if (is_aligner_0) {
      // auto vh_nearest = CGAL::nearest_neighbor(*VM_Voro<BaseV>::DT_, pit.get_vertex());
      // unsigned int idx1 = vh_nearest->info();
      // VM_Voro<BaseV>::v_arr_[idx0].collide(&VM_Voro<BaseV>::v_arr_[idx1], true, false);
    }
  }
}

template <class BaseV>
template <typename TRan>
void VM_Voro_AlignerDissenter<BaseV>::ini_rand(TRan &myran, double theta0) {
  double *x = new double[VM_Voro<BaseV>::N_];
  double *y = new double[VM_Voro<BaseV>::N_];
  double *vx = new double[VM_Voro<BaseV>::N_];
  double *vy = new double[VM_Voro<BaseV>::N_];

  double vx0 = cos(theta0);
  double vy0 = sin(theta0);
  for (int i = 0; i < VM_Voro<BaseV>::N_; i++) {
    x[i] = myran.doub() * VM_Voro<BaseV>::Lx_;
    y[i] = myran.doub() * VM_Voro<BaseV>::Ly_;
    if (i > n_dis_) {
      vx[i] = vx0;
      vy[i] = vy0;
    } else {
      double theta = myran.doub() * 2 * PI;
      vx[i] = cos(theta);
      vy[i] = sin(theta);
    }
  }
  VM_Voro<BaseV>::input_data(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
}

template <class BaseV>
template <typename TSnap>
void VM_Voro_AlignerDissenter<BaseV>::ini_from_snap(TSnap &snap_reader) {
  double *x = new double[VM_Voro<BaseV>::N_];
  double *y = new double[VM_Voro<BaseV>::N_];
  double *vx = new double[VM_Voro<BaseV>::N_];
  double *vy = new double[VM_Voro<BaseV>::N_];
  uint32_t *type_id = new uint32_t[VM_Voro<BaseV>::N_];
  snap_reader.read_last_frame(x, y, vx, vy, type_id);
  std::cout << "Load " << VM_Voro<BaseV>::N_ << " particles, consisting of " << n_dis_ << " dissenters" << std::endl;
  for(int i = 0; i < VM_Voro<BaseV>::N_; i++) {
    if (i < n_dis_) {
      if (type_id[i] != 1) {
        std::cout << "Error, typeid of particle " << i << " is " << type_id[i] << std::endl;
        exit(1);
      }
    } else {
      if (type_id[i] != 0) {
        std::cout << "Error, typeid of particle " << i << " is " << type_id[i] << std::endl;
        exit(1);
      }
    }
  }
  VM_Voro<BaseV>::input_data(x, y, vx, vy);
  delete []x;
  delete []y;
  delete []vx;
  delete []vy;
  delete []type_id;
}

template <class BaseV>
template <typename T>
void VM_Voro_AlignerDissenter<BaseV>::get_order_para(T &phi, T &theta) const {
    double vx = 0;
    double vy = 0;
    for (int i = n_dis_; i < VM_Voro<BaseV>::N_; i++)
    {
        vx += VM_Voro<BaseV>::v_arr_[i].vx;
        vy += VM_Voro<BaseV>::v_arr_[i].vy;
    }
    int n_aligner = VM_Voro<BaseV>::N_ - n_dis_;
    vx /= n_aligner;
    vy /= n_aligner;
    phi = sqrt(vx * vx + vy * vy);
    theta = atan2(vy, vx);
}

template <class BaseV>
template <typename TSnap>
void VM_Voro_AlignerDissenter<BaseV>::dump(int i_step, TSnap &snap_writer) {
  if (snap_writer.need_export(i_step)) {
    float* pos = new float[VM_Voro<BaseV>::N_*3];
    uint32_t *type_id = new uint32_t[VM_Voro<BaseV>::N_]{};
    for (int i = 0; i < n_dis_; i++) {
      type_id[i] = 1;
    }
    VM_Voro<BaseV>::get_pos_arr(pos);
    snap_writer.dump(i_step, VM_Voro<BaseV>::N_, pos, type_id);
    delete [] pos;
    delete [] type_id;
  }
}
*/