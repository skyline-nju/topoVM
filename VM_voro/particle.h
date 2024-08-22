#pragma once
#include <cmath>
#include "comn.h"

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

// velocity of particles that move in discrete time steps and suffer the 
// scalar noise.
struct V_scalar {
  V_scalar() {}

  V_scalar(double _vx, double _vy) :
    vx(_vx), vy(_vy), vx_next(_vx), vy_next(_vy) {
  }

  template <class T>
  void collide(T *other);

  template <class T>
  void collide(T *other, bool is_aligner_1, bool is_aligner_2);

  template <class T>
  void collide_asym(T *other, double dx, double dy, double alpha);

  template <typename TRan>
  void update_v(double eta, TRan &myran);

  template <typename T>
  void update_x(T &_x, T &_y, double v0, double Lx, double Ly);

  void get_v(double &VX, double &VY) const { VX = vx; VY = vy; }
  void get_theta(double &Theta) const { Theta = std::atan2(vy, vx); }

  double vx;
  double vy;
  double vx_next;
  double vy_next;
};

template <class T>
void V_scalar::collide(T *other) {
  vx_next += other->vx;
  vy_next += other->vy;
  other->vx_next += vx;
  other->vy_next += vy;
}

template <class T>
void V_scalar::collide(T *other, bool is_aligner_1, bool is_aligner_2) {
  if (is_aligner_1) {
    vx_next += other->vx;
    vy_next += other->vy;
  }
  if (is_aligner_2) {
    other->vx_next += vx;
    other->vy_next += vy;
  }
}

template <class T>
void V_scalar::collide_asym(T *other, double dx, double dy, double alpha) {
  double w1 = 0.5 * (1 + alpha * sgn(dx * vx + dy * vy));
  double w2 = 0.5 * (1 + alpha * sgn(-dx * other->vx - dy * other->vy));
  vx_next += w1 * other->vx;
  vy_next += w1 * other->vy;
  other->vx_next += w2 * vx;
  other->vy_next += w2 * vy;
}

template <typename TRan>
void V_scalar::update_v(double eta, TRan &myran) {
  double tmp = std::sqrt(vx_next * vx_next + vy_next * vy_next);
  double c1 = vx_next / tmp;
  double s1 = vy_next / tmp;
  double noise = (myran.doub() - 0.5) * eta * 2 * PI;

  double c2 = std::cos(noise);
  double s2 = std::sin(noise);
  vx = vx_next = c1 * c2 - s1 * s2;
  vy = vy_next = c1 * s2 + c2 * s1;
}

template <typename T>
void V_scalar::update_x(T & _x, T & _y, double v0, double Lx, double Ly) {
  _x += v0 * vx;
  _y += v0 * vy;
  if (_x < 0) {
    _x += Lx;
  } else if (_x >= Lx) {
    _x -= Lx;
  }
  if (_y < 0) {
    _y += Ly;
  } else if (_y >= Ly) {
    _y -= Ly;
  }
}