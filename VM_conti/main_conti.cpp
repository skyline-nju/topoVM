#include "vm_knn.h"
#include "io2D.h"

int main(int argc, char* argv[]) {
  // Set parameters
  double Lx = atof(argv[1]);
  double Ly = atof(argv[2]);
  double Dr = atof(argv[3]);
  double J0 = atof(argv[4]);
  int k = atoi(argv[5]); // number of neighbors for alignment
  double rho0 = 1.0;

  int n_step = atoi(argv[6]);
  int snap_dt = atoi(argv[7]);
  int seed = atoi(argv[8]);
  std::string ini_mode = argv[9];

  double h0 = 0.1;
  double v0 = 0.2;
  int N = int(Lx * Ly * rho0);
  int log_dt = 1000;
  Ranq2 myran(seed);

  double J0_over_k = J0 / k;
  VM_kNN birds(Lx, Ly, rho0, J0_over_k, Dr, h0, v0);

  // Set output
  char basename[255];
  char log_file[1024];
  char order_para_file[1024];
  char snap_file[1024];

  char prefix[200] = "./data";
  mkdir(prefix);
  char log_folder[255];
  char order_para_folder[255];
  snprintf(log_folder, 255, "%s/log", prefix);
  snprintf(order_para_folder, 255, "%s/op", prefix);
  mkdir(log_folder);
  mkdir(order_para_folder);


  snprintf(basename, 255, 
    "L%g_%g_k%d_D%.5f_J%.5f_r%.3f_h%.4f_s%d", 
    Lx, Ly, k, Dr, J0, rho0, h0, seed);
  snprintf(snap_file, 1024, "%s/%s.gsd", prefix, basename);

  int start = 0;
  double snap_log_sep = -1;
  Snap_GSD_2 gsd(snap_file, n_step, snap_dt, start, h0, snap_log_sep, Lx, Ly, ini_mode);

  snprintf(log_file, 1024, "%s/%s_t%08d.dat", log_folder, basename, start);
  snprintf(order_para_file, 1024, "%s/%s_t%08d.dat", order_para_folder, basename, start);

  LogExporter log(log_file, start, n_step, 1000, N);
  OrderParaExporter op(order_para_file, start, n_step, 100);

  if (ini_mode == "resume") {
    birds.ini_from_snap(gsd);
  } else if (ini_mode == "rand") {
    birds.ini_rand(myran);
  } else if (ini_mode == "ordered") {
    birds.ini_rand(myran, 0);
  } else {
    std::cout << "Error, ini_mode must be one of resume, rand, ordered and bi-ordered!" << std::endl;
    exit(1);
  }

  for (int i = 1; i <= n_step; i++) {
    birds.collide(k);
    birds.stream(myran);
    log.record(i);
    op.dump(i, birds);
    birds.dump(i, gsd);
  }
}

