// g++ main.cpp -std=c++17 -O3 -march=native -ffast-math -o solver_ref
// ./solver_ref reference 10 10000

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <iterator>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

struct CCSmatrix {
  std::vector<double> V;
  std::vector<size_t> IA;
  std::vector<size_t> JA;
  size_t size;
  size_t numEntries;

  void readMtxFile(std::ifstream &file) {
    int lineNum=0;
    std::string line;
    int currentCol = 0;
    while(getline(file, line)){
      //std::cout << line << std::endl;
      //std::cout << "lineNum: " << lineNum << std::endl;
      if (line[0] != '%') {
        std::istringstream iss(line);
        std::vector<double> nums;
        double value;
        while(iss>>value) {
          nums.push_back(value);
        }
        //if (lineNum < 50) {
          // for (auto const &n : nums) {
          //   std::cout << n << " ";
          // }
          // std::cout << std::endl;
        //}
        
        if (lineNum != 0) {
          IA[lineNum -1] = nums[0] - 1;
          if (currentCol != nums[1]) {
            JA[currentCol] = lineNum - 1;
            //if (lineNum < 50)
            //std::cout << "JA: " << JA[currentCol] << std::endl;
            currentCol++;
          }
          V[lineNum - 1] = nums[2];
        } else {
          size = nums[0];
          numEntries = nums[2];
          V.resize(numEntries);
          IA.resize(numEntries);
          JA.resize(size+1);
        }
        nums.clear();
        lineNum++;
      }
    }
  }
};

namespace program_options {

struct Options {
  std::string filename;
  size_t iters;
  void print() const {
    std::printf("filename: %s\n", filename.c_str());
    std::printf("iters: %zu\n", iters);
  }
};

auto parse(int argc, char *argv[]) {
  if (argc != 3)
    throw std::runtime_error("unexpected number of arguments");
  Options opts;
  opts.filename = argv[1];
  // TODO: Validate filename is available
  if (std::sscanf(argv[2], "%zu", &opts.iters) != 1 && opts.iters != 0)
    throw std::runtime_error("invalid parameter for iters");
  return opts;
}

} // namespace program_options

int main(int argc, char *argv[]) try {

  // parse args
  auto opts = program_options::parse(argc, argv);
  opts.print();
  auto mat = CCSmatrix();
  std::ifstream matrixFile;  
  matrixFile.open(opts.filename);
  mat.readMtxFile(matrixFile);

  // 2 norm
  // auto norm2 = [N = opts.N](const auto &vec) -> auto {
  //   double sum = 0.0;
  //   for (size_t j = 0; j < N; ++j)
  //     for (size_t i = 1; i < (N - 1); ++i)
  //       sum += vec[i + j * N] * vec[i + j * N];

  //   return std::sqrt(sum);
  // };

  // std::cout << "  norm2 = " << norm2(x2) << std::endl;

  return EXIT_SUCCESS;
} catch (std::exception &e) {
  std::cout << e.what() << std::endl;
  return EXIT_FAILURE;
}
