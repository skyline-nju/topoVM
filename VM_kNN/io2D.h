#pragma once
#include <vector>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "gsd.h"
#include "comn.h"

/**
 * @brief Basic class for exporting data.
 *
 * Define the timming to dump data.
 */
class ExporterBase {
public:
  ExporterBase(int n_step, int sep, int start) : n_step_(n_step), sep_(sep), start_(start) {}

  bool need_export(const int i_step);

  void set_log_scale_frames(double h, double log_sep = 0.1);

protected:
  int n_step_;    // total steps to run
  int sep_;
  int start_ = 0; // The first step 
  int my_rank_ = 0;
  int tot_proc_ = 1;
  std::vector<int> frames_;
  int cur_frame_ = 0;
  double log_sep_ = -1;
};

/**
 * @brief Exporter to output log
 *
 * Output the parameters after the initialization.
 * Output the beginning and endding time of the simulation.
 * Record time every certain time steps.
 */
class LogExporter : public ExporterBase {
public:
  LogExporter(const std::string& outfile, int start, int n_step, int sep,
    int np);

  ~LogExporter();

  void record(int i_step);

  std::ofstream fout;
private:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
  int n_par_;
  int step_count_ = 0;
};


class OrderParaExporter: public ExporterBase {
public:
  OrderParaExporter(const std::string& outfile, int start, int n_step, int sep)
    : ExporterBase(n_step, sep, start), fout_(outfile) {}

  ~OrderParaExporter() { fout_.close(); }

  template <typename TVM>
  void dump(int i_step, const TVM& birds);

private:
  std::ofstream fout_;
};

template <typename TVM>
void OrderParaExporter::dump(int i_step, const TVM &birds) {
  if (need_export(i_step)) {
    double phi, theta;
    birds.get_order_para(phi, theta);
    fout_ << i_step << "\t" << std::setprecision(8)
          << phi << "\t" << theta << "\n";
  }
}

class Snap_GSD_2 : public ExporterBase {
public:
  Snap_GSD_2(const std::string& filename,
             int n_step, int sep, int & start,
             double h, double log_sep,
             double Lx, double Ly,
             const std::string& open_flag);

  ~Snap_GSD_2();

  template <typename TPar>
  void get_data_from_par(const std::vector<TPar>& p_arr, float* pos) const;

  template <typename T>
  void get_x_y_theta(float* pos, T& x, T& y, T& theta) const;

  uint64_t get_time_step();

  int reset_start_time_step();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr, float* v = nullptr);

  void dump(int i_step, int n_par, const float* pos, const uint32_t* type_id = nullptr);

  template <typename TPar>
  void read(int i_frame, std::vector<TPar>& p_arr);

  template <typename T>
  void read(int i_frame, T *x, T *y, T *vx, T *vy, uint32_t *type_id = nullptr);

  template <typename TPar>
  void read_last_frame(std::vector<TPar>& p_arr);

  template <typename T>
  void read_last_frame(T* x, T* y, T* vx, T*vy, uint32_t *type_id = nullptr);

private:
  gsd_handle* handle_ = nullptr;
  double half_Lx_; 
  double half_Ly_;
  double Lx_;
  double Ly_;
};


template <typename TPar>
void Snap_GSD_2::get_data_from_par(const std::vector<TPar>& p_arr, float* pos) const {
  size_t n_par = p_arr.size();
  for (size_t j = 0; j < n_par; j++) {
    size_t j3 = j * 3;
    pos[j3    ] = p_arr[j].pos.x - half_Lx_;
    pos[j3 + 1] = p_arr[j].pos.y - half_Ly_;
    pos[j3 + 2] = p_arr[j].get_theta();
  }
}

template <typename TPar>
void Snap_GSD_2::dump(int i_step, const std::vector<TPar>& p_arr, float* v) {
  if (need_export(i_step)) {
    uint32_t n_par = p_arr.size();
    float* pos = new float[n_par * 3];
    get_data_from_par(p_arr, pos);
    //uint64_t step = get_time_step();
    uint64_t step = start_ + i_step;

    std::cout << "dump frame " << step << std::endl;
    gsd_write_chunk(handle_, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &step);
    gsd_write_chunk(handle_, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n_par);
    gsd_write_chunk(handle_, "particles/position", GSD_TYPE_FLOAT, n_par, 3, 0, pos);
    if (v != nullptr) {
      gsd_write_chunk(handle_, "particles/charge", GSD_TYPE_FLOAT, n_par, 1, 0, v);
    }
    gsd_end_frame(handle_);
    delete[] pos;
  }
}

template <typename TPar>
void Snap_GSD_2::read(int i_frame, std::vector<TPar>& p_arr) {
  uint32_t n_par;
  const gsd_index_entry* chunk = gsd_find_chunk(handle_, i_frame, "particles/N");
  gsd_read_chunk(handle_, &n_par, chunk);
  std::cout << "frame " << i_frame  <<": find " << n_par << " particles" << std::endl;

  float* pos = new float[n_par * 3];
  chunk = gsd_find_chunk(handle_, i_frame, "particles/position");
  gsd_read_chunk(handle_, pos, chunk);
  uint32_t* type_id = new uint32_t[n_par];

  p_arr.reserve(n_par);
  for (int j = 0; j < n_par; j++) {
    double x, y, theta;
    get_x_y_theta(pos + 3 * j, x, y, theta);

    TPar p;
    // p.pos = Vec_2<double>(x, y);
    // p.u = Vec_2<double>(cos(theta), sin(theta));
    // p_arr.push_back(p);
    p_arr.emplace_back(x, y, theta);
  }

  delete[] pos;
  delete[] type_id;
}


template <typename T>
void Snap_GSD_2::read(int i_frame, T *x, T *y, T *vx, T *vy, uint32_t *type_id) {
  uint32_t n_par;
  const gsd_index_entry* chunk = gsd_find_chunk(handle_, i_frame, "particles/N");
  gsd_read_chunk(handle_, &n_par, chunk);
  std::cout << "frame " << i_frame  <<": find " << n_par << " particles" << std::endl;
  float* pos = new float[n_par * 3];
  chunk = gsd_find_chunk(handle_, i_frame, "particles/position");
  gsd_read_chunk(handle_, pos, chunk);

  for (int j = 0; j < n_par; j++) {
    double theta;
    get_x_y_theta(pos + 3 * j, x[j], y[j], theta);
    vx[j] = cos(theta);
    vy[j] = sin(theta);
  }

  if (type_id) {
    chunk = gsd_find_chunk(handle_, i_frame, "particles/typeid");
    gsd_read_chunk(handle_, type_id, chunk);
  }
}


template <typename TPar>
void Snap_GSD_2::read_last_frame(std::vector<TPar>& p_arr) {
  int nframes = gsd_get_nframes(handle_);
  if (nframes < 1) {
    std::cout << "Error, nframes=" << nframes << std::endl;
    exit(1);
  } else {
    read(nframes-1, p_arr);
  }
}

template <typename T>
void Snap_GSD_2::read_last_frame(T *x, T *y, T *vx, T *vy, uint32_t *type_id) {
  int nframes = gsd_get_nframes(handle_);
  if (nframes < 1) {
    std::cout << "Error, nframes=" << nframes << std::endl;
    exit(1);
  } else {
    read(nframes-1, x, y, vx, vy, type_id);
  }
}

template <typename T>
void Snap_GSD_2::get_x_y_theta(float *pos, T &x, T &y, T &theta) const {
  x = pos[0] + half_Lx_;
  y = pos[1] + half_Ly_;
  theta = pos[2];
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

