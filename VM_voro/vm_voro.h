#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "particle.h"
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
 *      Class for one-species flocks with voronoi-type interaction
 * 
 * 
 ************************************************************************************/
template <class BaseV>
class VM_Voro: public VM {
public:
  VM_Voro(double Lx, double Ly, int N, double eta, double v0=0.5);

  ~VM_Voro(){ delete DT_; }

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

  void align();

  void align_nearest_neighbor();

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

  // array of vertex handle
  std::vector<Vertex_handle> vh_arr_;

  PDT *DT_;
};

template <class BaseV>
VM_Voro<BaseV>::VM_Voro(double Lx, double Ly, int N, double eta, double v0)
  : VM(Lx, Ly, N, eta, v0) {
  // initialize triangulation
  PDT::Iso_rectangle domain(0, 0, Lx_, Ly_);
  DT_ = new PDT(domain);
}

template<class BaseV>
template<typename T>
void VM_Voro<BaseV>::input_data(const T *x, const T *y, const T *vx, const T *vy) {
  v_arr_.reserve(N_);
  vh_arr_.reserve(N_);
  for (int i = 0; i < N_; i++) {
    v_arr_.emplace_back(vx[i], vy[i]);
    auto vh = DT_->insert(Point(x[i], y[i]));
    vh->info() = i;
    vh_arr_.push_back(vh);
  }
}

template <class BaseV>
template <typename T>
void VM_Voro<BaseV>::get_pos_arr(T *pos) const {
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
void VM_Voro<BaseV>::ini_rand(TRan &myran) {
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
void VM_Voro<BaseV>::ini_rand(TRan &myran, double theta0) {
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
void VM_Voro<BaseV>::ini_from_snap(TSnap &snap_reader) {
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
void VM_Voro<BaseV>::dump(int i_step, TSnap &snap_writer) {
  if (snap_writer.need_export(i_step)) {
    float* pos = new float[N_*3];
    get_pos_arr(pos);
    snap_writer.dum(i_step, N_, pos);
    delete [] pos;
  }
}

template <class BaseV>
void VM_Voro<BaseV>::align() {
  for (auto fit = DT_->periodic_triangles_begin(PDT::UNIQUE);
    fit != DT_->periodic_triangles_end(PDT::UNIQUE); ++fit) {
    unsigned int idx0 = fit.get_face()->vertex(0)->info();
    unsigned int idx1 = fit.get_face()->vertex(1)->info();
    unsigned int idx2 = fit.get_face()->vertex(2)->info();
    if (idx0 < idx1) {
      v_arr_[idx0].collide(&v_arr_[idx1]);
    }
    if (idx1 < idx2) {
      v_arr_[idx1].collide(&v_arr_[idx2]);
    }
    if (idx2 < idx0) {
      v_arr_[idx2].collide(&v_arr_[idx0]);
    }
  }
}

// TODO
template <class BaseV>
void VM_Voro<BaseV>::align_nearest_neighbor() {
  for (auto pit = DT_->periodic_points_begin(PDT::UNIQUE);
    pit != DT_->periodic_points_end(PDT::UNIQUE); ++pit) {
    unsigned int idx0 = pit.get_vertex()->info();
    // auto vh_nearest = CGAL::nearest_neighbor(*DT_, pit.get_vertex());
    // unsigned int idx1 = vh_nearest->info();
    // v_arr_[idx0].collide(&v_arr_[idx1], true, false);
  }
}

template <class BaseV>
template <typename TRan>
void VM_Voro<BaseV>::stream(TRan &myran) {
  for (int i = 0; i < N_; i++) {
    v_arr_[i].update_v(eta_, myran);
    double x, y;
    get_x(i, x, y);
    v_arr_[i].update_x(x, y, v0_, Lx_, Ly_);
    Point p_new(x, y);
    auto vh_old = vh_arr_[i];
    vh_arr_[i] = DT_->insert(p_new, vh_old->face());
    vh_arr_[i]->info() = vh_old->info();
    DT_->remove(vh_old);
  }
}

template<class BaseV>
void VM_Voro<BaseV>::get_x(int i, double & x, double & y) const {
  auto p = vh_arr_[i]->point();
  x = p.x();
  y = p.y();
}

template <class BaseV>
template <typename T>
void VM_Voro<BaseV>::get_order_para(T &phi, T &theta) const {
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