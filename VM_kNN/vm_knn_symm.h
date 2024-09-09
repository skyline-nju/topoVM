#pragma once
#include "vm_knn.h"
#include <CGAL/Orthogonal_incremental_neighbor_search.h>

typedef CGAL::Orthogonal_incremental_neighbor_search<Traits> NN_incremental_search;
typedef NN_incremental_search::iterator NN_iterator;
typedef NN_incremental_search::Tree Tree_incremental;


// A functor that returns true, iff a point is not one of the k/2 front NNs or
// k/2 back NNs of particle i.
template <typename BaseV>
class X_not_half_k_NN {
public:
  X_not_half_k_NN(const Point_2& p0, double vx0, double vy0, int half_k)
    : p0_(p0), vx0_(vx0), vy0_(vy0), half_k_(half_k) {}

  bool operator()(const NN_iterator& it);

  bool stop() {return n_front_ >= half_k_ && n_back_ >= half_k_; }

  void show_nn();

private:
  int n_front_ = 0;
  int n_back_ = 0;

  Point_2 p0_;
  double vx0_;
  double vy0_;
  int half_k_;
};

template <typename BaseV>
bool X_not_half_k_NN<BaseV>::operator()(const NN_iterator &it) {
  bool res = true;
  const Point_2 &p1 =  boost::get<0>(it->first);

  if (p0_ != p1) {
    double dx = p1.x() - p0_.x();
    double dy = p1.y() - p0_.y(); 
    
    double cross_dot = vx0_ * dx + vy0_ * dy;
    if (cross_dot < 0) {
      if (n_back_ < half_k_) {
        res = false;
      }
      n_back_++;
      // std::cout << "n_back_=" << n_back_ << std::endl;
    } else {
      if (n_front_ < half_k_) {
        res = false;
      }
      n_front_++;
      // std::cout << "n_front_=" << n_front_ << std::endl;
    }
    // show_nn();
  }
  
  return res; 
}

template <typename BaseV>
void X_not_half_k_NN<BaseV>::show_nn() {
  std::cout << "front neighbors: " << n_front_ << ", back neighbors: " << n_back_ << std::endl;
}

template <typename BaseV>
class X_not_front_nn {
public:
  X_not_front_nn(const Point_2& p0, double vx0, double vy0)
    : p0_(p0), vx0_(vx0), vy0_(vy0) {}

  bool operator()(const NN_iterator& it);

private:
  Point_2 p0_;
  double vx0_;
  double vy0_;
};

template <typename BaseV>
bool X_not_front_nn<BaseV>::operator()(const NN_iterator& it) {
  bool res = true;
  const Point_2 &p1 =  boost::get<0>(it->first);
  if (p0_ != p1) {
    double dx = p1.x() - p0_.x();
    double dy = p1.y() - p0_.y();
    res = vx0_ * dx + vy0_ * dy < 0;
  }
  return res;
}


template <typename BaseV>
class X_not_back_nn {
public:
  X_not_back_nn(const Point_2& p0, double vx0, double vy0)
    : p0_(p0), vx0_(vx0), vy0_(vy0) {}

  bool operator()(const NN_iterator& it);

private:
  Point_2 p0_;
  double vx0_;
  double vy0_;
};

template <typename BaseV>
bool X_not_back_nn<BaseV>::operator()(const NN_iterator &it) {
  bool res = true;
  const Point_2 &p1 =  boost::get<0>(it->first);
  if (p0_ != p1) {
    double dx = p1.x() - p0_.x();
    double dy = p1.y() - p0_.y();
    res = vx0_ * dx + vy0_ * dy > 0;
  }
  return res;
}



/***********************************************************************************
 *   Class for one-species flocks with kNN interaction with fore-aft symmetry,
 *   such that one particle looks for its n/2 front neighbors and n/2 back neighbors.
 ***********************************************************************************/
template <class BaseV>
class VM_kNN_symm: public VM_kNN<BaseV> {
public:
  typedef CGAL::Filter_iterator<NN_iterator, X_not_half_k_NN<BaseV>> NN_filter_iterator;
  typedef CGAL::Filter_iterator<NN_iterator, X_not_front_nn<BaseV>> front_NN_iterator;
  typedef CGAL::Filter_iterator<NN_iterator, X_not_back_nn<BaseV>> back_NN_iterator;

  VM_kNN_symm(double Lx, double Ly, int N, double eta, double v0=0.5)
    : VM_kNN<BaseV>(Lx, Ly, N, eta, v0) {}

  void align_base(int k, double *pt_max_dis2=nullptr);

  void align(int k);

  void align_padded(int k);
};

template <class BaseV>
void VM_kNN_symm<BaseV>::align_base(int k, double *pt_max_dis2) {
  Tree tree(
    boost::make_zip_iterator(boost::make_tuple( VM_kNN<BaseV>::pos_arr_.begin(),VM_kNN<BaseV>::idx_arr_.begin())),
    boost::make_zip_iterator(boost::make_tuple( VM_kNN<BaseV>::pos_arr_.end(),VM_kNN<BaseV>::idx_arr_.end())));
  
  for (int i = 0; i < VM_kNN<BaseV>::N_; i++) {
    double vx0, vy0;
    VM_kNN<BaseV>::get_v(i, vx0, vy0);
    Point_2 p0 = VM_kNN<BaseV>::pos_arr_[i];
    NN_incremental_search NN(tree, p0);

    // the first k/2 front neighbors of particle i
    { 
      X_not_front_nn<BaseV> filter(p0, vx0, vy0);
      front_NN_iterator it(NN.end(), filter, NN.begin());
      front_NN_iterator end(NN.end(), filter);
      int count = 0;
      while (count < k/2 && it != end) {
        int j = boost::get<1>(it->first);
        VM_kNN<BaseV>::v_arr_[i].collide_one_side(VM_kNN<BaseV>::v_arr_[j]);
        if (pt_max_dis2 && *pt_max_dis2 < it->second) {
          *pt_max_dis2 = it->second;
        }
        ++it;
        count++;
        // std::cout << "j=" << j << ", d2=" << it->second << std::endl;
      }
    }

    // the first k/2 back neighbors of particle i
    { 
      X_not_back_nn<BaseV> filter(p0, vx0, vy0);
      back_NN_iterator it(NN.end(), filter, NN.begin());
      back_NN_iterator end(NN.end(), filter);
      int count = 0;
      while (count < k/2 && it != end) {
        int j = boost::get<1>(it->first);
        VM_kNN<BaseV>::v_arr_[i].collide_one_side(VM_kNN<BaseV>::v_arr_[j]);
        if (pt_max_dis2 && *pt_max_dis2 < it->second) {
          *pt_max_dis2 = it->second;
        }
        ++it;
        count++;
        // std::cout << "j=" << j << ", d2=" << it->second << std::endl;
      }
    }
    // exit(1);
  }
}

template <class BaseV>
void VM_kNN_symm<BaseV>::align(int k) {
  VM_kNN<BaseV>::add_padded_particles();
  align_base(k);
  VM_kNN<BaseV>::del_padded_particles();
}

template <class BaseV>
void VM_kNN_symm<BaseV>::align_padded(int k) {
  if (VM_kNN<BaseV>::padded_dx2_ == 0 || VM_kNN<BaseV>::padded_dx2_ >= VM_kNN<BaseV>::Lx_ * VM_kNN<BaseV>::Lx_ || VM_kNN<BaseV>::padded_dx2_ >= VM_kNN<BaseV>::Ly_ * VM_kNN<BaseV>::Ly_) {
    VM_kNN<BaseV>::add_padded_particles();
    double dx2 = 0;
    align_base(k, &dx2);
    VM_kNN<BaseV>::padded_dx2_ = dx2;
    VM_kNN<BaseV>::del_padded_particles();
  } else {
    double dx = sqrt(VM_kNN<BaseV>::padded_dx2_) + 1;
    double dy = dx;
    VM_kNN<BaseV>::add_padded_particles(dx, dy);
    double dx2 = 0;
    align_base(k, &dx2);
    VM_kNN<BaseV>::padded_dx2_ = dx2;
    VM_kNN<BaseV>::del_padded_particles();
  }
}
