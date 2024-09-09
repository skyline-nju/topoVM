#pragma once
#include "vm_knn.h"
#include <CGAL/Orthogonal_incremental_neighbor_search.h>

typedef CGAL::Orthogonal_incremental_neighbor_search<Traits> NN_incremental_search;
typedef NN_incremental_search::iterator NN_iterator;
typedef NN_incremental_search::Tree Tree_incremental;

/***********************************************************************************
 *   Class for one-species flocks with kNN interaction with fore-aft symmetry,
 *   such that one particle looks for its n/2 front neighbors and n/2 back neighbors.
 ***********************************************************************************/
template <class BaseV>
class VM_kNN_symm: public VM_kNN<BaseV> {
public:
  VM_kNN_symm(double Lx, double Ly, int N, double eta, double v0=0.5)
    : VM_kNN<BaseV>(Lx, Ly, N, eta, v0) {}

  void align_base(int k, double *pt_max_dis2=nullptr);

  void align(int k);

  void align_padded(int k);
};

template <class BaseV>
void VM_kNN_symm<BaseV>::align_base(int k, double *pt_max_dis2) {
  Tree_incremental tree(
    boost::make_zip_iterator(boost::make_tuple( VM_kNN<BaseV>::pos_arr_.begin(),VM_kNN<BaseV>::idx_arr_.begin())),
    boost::make_zip_iterator(boost::make_tuple( VM_kNN<BaseV>::pos_arr_.end(),VM_kNN<BaseV>::idx_arr_.end())));
  
  for (int i = 0; i < VM_kNN<BaseV>::N_; i++) {
    double vx0, vy0;
    VM_kNN<BaseV>::get_v(i, vx0, vy0);
    Point_2 p0 = VM_kNN<BaseV>::pos_arr_[i];
    NN_incremental_search NN(tree, p0);
    int half_k = k / 2;
    int front_nn = 0;
    int back_nn = 0;

    for (auto it=NN.begin(); it!=NN.end();++it) {
      int j = boost::get<1>(it->first);
      if (j != i) {
        const Point_2 &p1 =  boost::get<0>(it->first);
        double dx = p1.x() - p0.x();
        double dy = p1.y() - p0.y();
        if (vx0 * dx + vy0 * dy < 0) {
          if (back_nn < half_k) {
            VM_kNN<BaseV>::v_arr_[i].collide_one_side(VM_kNN<BaseV>::v_arr_[j]);
            // std::cout << "back NN with d2 " << it->second << std::endl;
          }
          back_nn++;
        } else {
          if (front_nn < half_k) {
            VM_kNN<BaseV>::v_arr_[i].collide_one_side(VM_kNN<BaseV>::v_arr_[j]);
            // std::cout << "front NN with d2 " << it->second << std::endl;
          }
          front_nn++;
        }
        if (pt_max_dis2 && *pt_max_dis2 < it->second) {
          *pt_max_dis2 = it->second;
        }
      }
      if (back_nn >= half_k && front_nn >= half_k) {
        break;
      } 
    }
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
