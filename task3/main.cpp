// g++ main.cpp -std=c++17 -O3 -march=native -ffast-math -o solver

#include <cmath>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <ctime>
#include <random>
#include <omp.h>
#include <chrono>
#include <time.h>
//#include <functional>

const double PI = 3.141592653589793238463;

// Structure to hold seed with proper cache line alignment
struct alignas(64) ThreadSeed {
  unsigned int seed;
  char padding[60];
};

namespace program_options {

struct Options {
  std::string func;
  double x_min;
  double x_max;
  int samples;

  std::unordered_set<std::string> valid_funcs = {"SINX", "COS2XINV", "X4M5"};

  void print() const {
    std::printf("func: %s\n", func.c_str());
    std::printf("x_min: %lf\n", x_min);
    std::printf("x_max: %lf\n", x_max);;
    std::printf("samples: %d\n", samples);
  }
};

auto parse(int argc, char *argv[]) {
  if (argc != 5)
    throw std::runtime_error("unexpected number of arguments");
  Options opts;
  opts.func = argv[1];
  if (opts.valid_funcs.find(opts.func) == opts.valid_funcs.end()) {
    std::stringstream error_message;
    error_message << "invalid parameter for func, valid values are: ";
    for (const auto& func : opts.valid_funcs) {
      error_message << func << " ";
    }
    throw std::runtime_error(error_message.str());
  }
  if (std::sscanf(argv[2], "%lf", &opts.x_min) != 1)
    throw std::runtime_error("invalid parameter for x_min");
  if (std::sscanf(argv[3], "%lf", &opts.x_max) != 1)
    throw std::runtime_error("invalid parameter for x_max");
  double temp = 0;

  if (std::sscanf(argv[4], "%lf", &temp) != 1 && temp <= 0)
    throw std::runtime_error("invalid parameter for samples");
  else
    opts.samples = static_cast<int>(temp);
  return opts;
}

} // namespace program_options

struct Integrator {
  double x_min, x_max;
  int samples;
  std::string func;
  std::mt19937_64 mt;
  //std::function<double(double)> selectedFunc;
  double (*selectedFunc)(double) = nullptr;

  static double sinx(double x){
    return sin(x);
  }

  static double cos2xinv(double x){
    double result = cos(1/x)*cos(1/x);
    return result;
  }

  static double x4m5(double x){
    double result = 5*x*x*x*x;
    return result;
  }

  void setup(program_options::Options& options){
    func = options.func;
    if(func == "SINX") {
      selectedFunc = sinx;
    }
    else if(func == "COS2XINV")
      selectedFunc = cos2xinv;
    else if(func == "X4M5")
      selectedFunc = x4m5;
      
    x_min = options.x_min;
    x_max = options.x_max;
    samples = options.samples;
  }

  

  double scaler(double x){
    double result = x_min + (x_max-x_min)/mt.max()*x;
    return result;
  }

  double integrate(){
    mt.seed(time(nullptr));
    double sum = 0;
  
    #pragma omp parallel shared (sum)
    {
    std::cout << "thread id:" << omp_get_thread_num() << std::endl;
    #pragma omp parallel for reduction(+:sum)
      for(int i = 0; i < samples; i++){
        double x = scaler(mt());
        sum += selectedFunc(x); 
      }
    }
    sum = sum/samples*(x_max-x_min);
    return sum;
  }

  
};



int main(int argc, char *argv[]) try {
  //int numThreads = omp_get_max_threads();
  //omp_set_num_threads(numThreads);
  std::cout << "max threads: " << omp_get_max_threads()<< std::endl;
  std::cout << "num threads: " << omp_get_num_threads()<< std::endl;
  bool debug = false;
  auto opts = program_options::parse(argc, argv);
  opts.print();

  // Get maximum number of available threads
  int numThreads = omp_get_max_threads();

  // Set the number of threads
  omp_set_num_threads(numThreads);

  // Initialization of seed "bank" (one seed per thread)
  std::vector<ThreadSeed> seedBank(numThreads);
  
  // Shared RNG to generate the seeds of each thread
  // std::random_device rd;
  // unsigned int seed = rd();
  unsigned int seed = 12347308;
  std::mt19937_64 sharedGen(12347308);
  std::uniform_int_distribution<unsigned int> dis;
  for (int i = 0; i < numThreads; i++) {
    seedBank[i].seed = dis(sharedGen);
    std::cout << seedBank[i].seed << std::endl;
  }



 
  // auto integrator_m = Integrator();
  // integrator_m.setup(opts);
  // double time_s = omp_get_wtime();
  // struct timeval start, end;
  // gettimeofday(&start, NULL);
  // std::cout<< integrator_m.integrate() << std::endl;
  // gettimeofday(&end, NULL);
  // delta = ((end.tv_sec-start.tv_sec)*1000000u + end.tv_sec-start.tv_sec)/1e6 << std::endl;
  // double time_e = omp_get_wtime();
  // std::cout << "time:" << (time_e-time_s)*1e3 << std::endl;
    

  return EXIT_SUCCESS;
} catch (std::exception &e) {
  std::cout << e.what() << std::endl;
  return EXIT_FAILURE;
}
