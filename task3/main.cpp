#include <cmath>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <chrono>
#include <random>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <chrono>

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

double sinx(double x) {
  return sin(x);
}

double cos2xinv(double x) {
  double cos1qx = cos(1/x);
  return cos1qx*cos1qx;
}

double x4m5(double x) {
  return 5*x*x*x*x;
}

int main(int argc, char *argv[]) try {
  bool debug = false;
  auto opts = program_options::parse(argc, argv);
  if(debug) opts.print();
  
  // Set selected function from opts
  double (*selectedFunc)(double) = nullptr;
  if (opts.func == "SINX")
    selectedFunc = sinx;
  else if (opts.func == "COS2XINV")
    selectedFunc = cos2xinv;
  else
    selectedFunc = x4m5;

  // Get maximum number of available threads
  int num_threads = 1;
  #ifdef _OPENMP
    num_threads = omp_get_max_threads();
    // Set the number of threads
    omp_set_num_threads(num_threads);
  #endif

  // Initialization of seed "bank" (one seed per thread)
  std::vector<ThreadSeed> seed_bank(num_threads);
  
  // SHARED RNG to generate the seeds of each thread
  // std::random_device rd;
  // unsigned int seed = rd();
  unsigned int seed = 12347308; // fixed seed
  std::mt19937_64 shared_gen(seed);
  std::uniform_int_distribution<unsigned int> int_dis;
  for (int i = 0; i < num_threads; i++) {
    seed_bank[i].seed = int_dis(shared_gen);
  }

  double sum = 0.0;
  
  // Start time measurement
  
  #ifdef _OPENMP
    double start_omp_time = omp_get_wtime();
  #else
    auto start_time = std::chrono::high_resolution_clock::now();
  #endif
  
  #pragma omp parallel reduction(+:sum)
  {
    int thread_num = 0;
    #ifdef _OPENMP
      thread_num = omp_get_thread_num();
    #endif

    std::mt19937_64 local_gen(seed_bank[thread_num].seed);
    std::uniform_real_distribution<double> real_dist(opts.x_min, opts.x_max);
    
    int count = 0; // Used calculate the integral approximation inside each thread

    #pragma omp for
    for(int i = 0; i < opts.samples; i++) {
      double x = real_dist(local_gen);
      double y = selectedFunc(x);
      sum += y;
      count++;
    }
    double thread_integral = sum/count*(opts.x_max-opts.x_min);
    printf("(%d)\t%f\n", thread_num, thread_integral);
  }
  if(debug) printf("total sum = %f\n", sum);
  double mean = sum/opts.samples;
  if(debug) printf("mean = %f\n", mean);
  double final_integral = mean*(opts.x_max-opts.x_min);
  printf("Integral:\t%f\n", final_integral);

  // End time measurement
  #ifdef _OPENMP
    double end_omp_time = omp_get_wtime();
    double elapsed_omp_time = end_omp_time - start_omp_time;
  #else
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time;
  #endif

  #ifdef _OPENMP
    printf("Runtime:\t%f\n",elapsed_omp_time);
  #else
    printf("Runtime:\t%f\n",elapsed_time.count());
  #endif

  printf("Threads:\t%d\n", num_threads);
  printf("Samples:\t%d\n", opts.samples);

  return EXIT_SUCCESS;
} catch (std::exception &e) {
  std::cout << e.what() << std::endl;
  return EXIT_FAILURE;
}
