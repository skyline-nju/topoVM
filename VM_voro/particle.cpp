#include "particle.h"

double V_conti::h;
double V_conti::J0_h;
double V_conti::sqrt_24_Dr_h;

void V_scalar::update_x(double & _x, double & _y, double v0, double Lx, double Ly) {
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

void V_conti::update_x(double & _x, double & _y, double v0, double Lx, double Ly) {
  double vx = v0 * std::cos(theta);
  double vy = v0 * std::sin(theta);
  _x += vx * h;
  _y += vy * h;
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

void V_conti::set_para(double h0, double Dr, double J0) {
  h = h0;
  sqrt_24_Dr_h = sqrt(24. * Dr * h);
  J0_h = h * J0;
}