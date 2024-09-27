#include "config.h"
#ifdef DISCRET_T
#include "vm_voro.h"
#include "io2D.h"

int main(int argc, char* argv[]) {
  // Set parameters
  double Lx = atof(argv[1]);
  double eta = atof(argv[2]);
  double frac_dis = atof(argv[3]);
  int n_step = atoi(argv[4]);
  int snap_dt = atoi(argv[5]);
  int seed = atoi(argv[6]);
  std::string ini_mode = argv[7];
  double Ly = Lx;
  double rho0 = 1;
  double v0 = 0.5;
  int N = int(Lx * Ly * rho0);
  int log_dt = 1000;
  int n_dis = int(round(N * frac_dis));
  Ranq2 myran(seed);

  VM_Voro_AlignerDissenter<V_scalar> birds(Lx, Ly, N, eta, v0, n_dis);

  // Set output
  char basename[255];
  char log_file[255];
  char order_para_file[255];
  char snap_file[255];

  char prefix[255] = "data";
  mkdir(prefix);
  char log_folder[255];
  char order_para_folder[255];
  snprintf(log_folder, 255, "%s/log", prefix);
  snprintf(order_para_folder, 255, "%s/op", prefix);
  mkdir(log_folder);
  mkdir(order_para_folder);

  snprintf(basename, 255, "L%g_%g_d%.4f_e%.3f_r%g_s%d", Lx, Ly, frac_dis, eta, rho0, seed);
  snprintf(snap_file, 255, "%s/%s.gsd", prefix, basename);

  int start = 0;
  double h = 1;
  double snap_log_sep = -1;
  Snap_GSD_2 gsd(snap_file, n_step, snap_dt, start, h, snap_log_sep, Lx, Ly, ini_mode);

  snprintf(log_file, 255, "%s/%s_t%08d.dat", log_folder, basename, start);
  snprintf(order_para_file, 255, "%s/%s_t%08d.dat", order_para_folder, basename, start);

  LogExporter log(log_file, start, n_step, 1000, N);
  OrderParaExporter op(order_para_file, start, n_step, 100);

  if (ini_mode == "resume") {
    birds.ini_from_snap(gsd);
  } else if (ini_mode == "rand") {
    birds.ini_rand(myran);
  } else if (ini_mode == "ordered") {
    birds.ini_rand(myran, 0);
  } else {
    std::cout << "Error, ini_mode must be one of resume, rand and ordered!" << std::endl;
    exit(1);
  }

  for (int i = 1; i <= n_step; i++) {
    birds.align();
    birds.stream(myran);
    log.record(i);
    op.dump(i, birds);
    birds.dump(i, gsd);
  }
}
#endif
