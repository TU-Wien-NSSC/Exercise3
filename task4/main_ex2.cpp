// g++ main.cpp -std=c++17 -O3 -march=native -ffast-math -o solver
// ./solver without_source 10 10000 0.0 1.0
// ./solver with_source 50 100000 0.0 0.0 0.5 0.5 0.1

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <chrono>

const double PI = 3.141592653589793238463;

namespace program_options {

struct Options {
  std::string name;
  size_t N;
  size_t iters;
  double fix_west;
  double fix_east;
  bool has_source;
  double source_x;
  double source_y;
  double source_sigma;
  void print() const {
    std::printf("name: %s\n", name.c_str());
    std::printf("N: %zu\n", N);
    std::printf("iters: %zu\n", iters);
    std::printf("fix_west: %lf\n", fix_west);
    std::printf("fix_east: %lf\n", fix_east);
    std::printf("has_source: %s\n", has_source ? "true" : "false");
    std::printf("source_x: %lf\n", source_x);
    std::printf("source_y: %lf\n", source_y);
    std::printf("source_sigma: %lf\n", source_sigma);
  }
};

auto parse(int argc, char *argv[]) {
  if (argc != 9 && argc != 6)
    throw std::runtime_error("unexpected number of arguments");
  Options opts;
  opts.name = argv[1];
  if (std::sscanf(argv[2], "%zu", &opts.N) != 1 && opts.N >= 2)
    throw std::runtime_error("invalid parameter for N");
  if (std::sscanf(argv[3], "%zu", &opts.iters) != 1 && opts.iters != 0)
    throw std::runtime_error("invalid parameter for iters");
  if (std::sscanf(argv[4], "%lf", &opts.fix_west) != 1)
    throw std::runtime_error("invalid value for fix_west");
  if (std::sscanf(argv[5], "%lf", &opts.fix_east) != 1)
    throw std::runtime_error("invalid value for fix_east");
  if (argc == 6) {
    opts.has_source = false;
    opts.source_x = NAN;
    opts.source_y = NAN;
    opts.source_sigma = NAN;
    return opts;
  }
  if (std::sscanf(argv[6], "%lf", &opts.source_x) != 1)
    throw std::runtime_error("invalid value for source_x");
  if (std::sscanf(argv[7], "%lf", &opts.source_y) != 1)
    throw std::runtime_error("invalid value for source_y");
  if (std::sscanf(argv[8], "%lf", &opts.source_sigma) != 1)
    throw std::runtime_error("invalid value for source_sigma");
  opts.has_source = true;
  return opts;
}

} // namespace program_options

namespace solver {
struct System {
  // simulation name
  std::string name; 

  // (square) domain side length
  double L;

  // mandatory arguments
  int N;
  int iters;
  double uW;
  double uE;
  bool withSource;

  // optional arguments
  double xSrc;
  double ySrc;
  double sigmaSrc;

  // mesh size
  double h;

  // system's number of equations will be (N*2)^2 when variables are initialized
  int numEqs;

  // data structures for the system of equations
  std::vector<double> uh; // solution
  std::vector<double> bh; // right hand side
  std::vector<double> res; // residual

  // runtime
  double runtime = 0.0;
  
  void setup(program_options::Options& options);
  void build();
  void solve();
  void resid();
  double norm2();
  double normInf();
  void write();
  void print() {
    std::cout << "bh = " << std::endl;
    for (int k = 0; k < numEqs; k++) {
      std::cout << bh[k] << std::endl;
    }
    std::cout << "uh = " << std::endl;
    for (int k = 0; k < numEqs; k++) {
      std::cout << uh[k] << std::endl;
    }
    std::cout << "res = " << std::endl;
    for (int k = 0; k < numEqs; k++) {
      std::cout << res[k] << std::endl;
    }
  }
};

} //namespace solver

int main(int argc, char *argv[]) try {
  bool debug = false;
  auto opts = program_options::parse(argc, argv);
  // opts.print();

  auto system = new solver::System();
  system->setup(opts);
  system->build();
  system->solve();
  system->resid();
  if (debug) system->print();
  std::cout << system->norm2() << std::endl;
  std::cout << system->normInf() << std::endl;
  std::cout << system->runtime << std::endl;
  system->write();

  return EXIT_SUCCESS;
} catch (std::exception &e) {
  std::cout << e.what() << std::endl;
  return EXIT_FAILURE;
}

void solver::System::setup(program_options::Options& options) {
  // set values to options
  name  = options.name;
  N     = options.N;
  iters = options.iters;
  uW    = options.fix_west;
  uE    = options.fix_east;
  
  withSource = false;
  if (options.has_source) {
    withSource = true;
    xSrc       = options.source_x;
    ySrc       = options.source_y;
    sigmaSrc   = options.source_sigma;
  }
}

void solver::System::build() {
  // domain side length
  L = 1.0;

  // mesh size
  h = L/(N-1);

  // number of equations
  numEqs = (N-2)*(N-2);
  
  // initialize vectors uh, bh, res
  uh = std::vector<double>(numEqs,0.0);
  bh = std::vector<double>(numEqs,0.0);
  res = std::vector<double>(numEqs,0.0);

  double x,y;
  int i,j,k;

  if (withSource){
    double twoSigmaSquared = 2*sigmaSrc*sigmaSrc;
    // We rescale the system by h^2, i.e. multiply Au=b by h^2, to avoid that operation in the solver (when calculating the residual, we'll need to rescale back)
    double constantTerm = h*h/(PI*twoSigmaSquared);
    for(k = 1; k <= N-2; k++) {
      for (j = 1; j <= N-2; j++) {
        i = j + (k-1)*(N-2) - 1;
        x = j*h;
        y = k*h;
        bh[i] = -constantTerm*exp(-((x-xSrc)*(x-xSrc)+(y-ySrc)*(y-ySrc))/(twoSigmaSquared));
      }
    }
  }
}


void solver::System::solve() {
  auto uhOld = std::vector<double>(numEqs,0.0);
  auto uhNew = std::vector<double>(numEqs,0.0);

  double total_runtime = 0.0;
  auto start = std::chrono::steady_clock::now(); // Start time

  for (int l = 0; l < iters; l++) {
    for (int k = 1; k <= (N-2); k++) {
      for (int j = 1; j <= (N-2); j++) {
        int i  = j + (k-1)*(N-2) - 1 ;

        double sum = 0.0;

        double left  = j > 1     ? uhOld[i-1]     : uW;
        double right = j < (N-2) ? uhOld[i+1]     : uE;
        double down  = k > 1     ? uhOld[i-(N-2)] : 0.0;
        double up    = k < (N-2) ? uhOld[i+(N-2)] : 0.0;

        sum = down + left + right + up;
        
        double Aii = (k == 1 || k == (N-2)) ? -3.0 : -4.0;
        
        uhNew[i] = (bh[i] - sum)/Aii;
      }
    }
    if (l < iters - 1){
      std::swap(uhOld, uhNew);
    }
  }
  auto end = std::chrono::steady_clock::now(); // End time
  total_runtime += std::chrono::duration<double>(end - start).count(); // Accumulate runtime
  
  std::swap(uh,uhNew);
  runtime = total_runtime;
}

void solver::System::resid() {
  for (int k = 1; k <= (N-2); k++) {
    for (int j = 1; j <= (N-2); j++) {
      int i  = j + (k-1)*(N-2) - 1 ;
      double sum = 0.0;

      double left   = j > 1     ? uh[i-1]     : uW;
      double right  = j < (N-2) ? uh[i+1]     : uE;
      double center = (k == 1 || k == (N-2)) ? -3.0*uh[i] : -4.0*uh[i];
      double down   = k > 1     ? uh[i-(N-2)] : 0.0;
      double up     = k < (N-2) ? uh[i+(N-2)] : 0.0;

      sum = down + left + center + right + up;
      
      res[i] = (sum - bh[i])/(h*h); // Rescale by h^2 to obtain the residual of the original system Au=b (see: build() function)
    }
  }
}

double solver::System::norm2() {
  double sum = 0.0;
  for (int i = 0; i < numEqs; i++) {
    sum += res[i]*res[i];
  }
  return sqrt(sum);
}

double solver::System::normInf() {
  double norm = 0.0;
  for (int i = 0; i < numEqs; i++) {
    if (norm < fabs(res[i])) {
      norm = fabs(res[i]);
    }
  }
  return norm;
}

void solver::System::write() {
    std::ofstream csv;
    csv.open(name + ".csv");
    for (int j = 0; j < N-2; ++j) {
      for (int i = 0; i < N-2; ++i) {
        csv << uh[i + j * (N-2)];
        if (i != N-3) {
          csv << " ";
        }
      }
      csv << "\n";
    }
    csv.close();
}