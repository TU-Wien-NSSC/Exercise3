// g++ main.cpp -std=c++17 -O3 -march=native -ffast-math -o solver

#include <cmath>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>

const double PI = 3.141592653589793238463;

namespace program_options {

struct Options {
  std::string func;
  double x_min;
  double x_max;
  size_t samples;

  std::unordered_set<std::string> valid_funcs = {"SINX", "COS2XINV", "X4M5"};

  void print() const {
    std::printf("func: %s\n", func.c_str());
    std::printf("x_min: %lf\n", x_min);
    std::printf("x_max: %lf\n", x_max);;
    std::printf("samples: %zu\n", samples);
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
  if (std::sscanf(argv[3], "%lf", &opts.x_min) != 1)
    throw std::runtime_error("invalid parameter for x_max");
  if (std::sscanf(argv[4], "%zu", &opts.samples) != 1 && opts.samples != 0)
    throw std::runtime_error("invalid parameter for samples");
  return opts;
}

} // namespace program_options

int main(int argc, char *argv[]) try {
  bool debug = false;
  auto opts = program_options::parse(argc, argv);
  opts.print();

  return EXIT_SUCCESS;
} catch (std::exception &e) {
  std::cout << e.what() << std::endl;
  return EXIT_FAILURE;
}
