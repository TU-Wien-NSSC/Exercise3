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
#include <numeric>

struct CCSmatrix {
  std::vector<double> V;
  std::vector<size_t> IA;
  std::vector<size_t> JA;
  size_t size;
  size_t numEntries;

  void readMtxFile(std::ifstream &file) { 
    int lineNum = 0;
    std::string line;
    int currentCol = 0;
    
    while (getline(file, line)) {
      if (line[0] != '%') { 
        std::istringstream iss(line);
        std::vector<double> nums;
        double value;
        while (iss >> value) {
          nums.push_back(value);
        }

        if (lineNum != 0) {
            size_t currentRow = (nums[0] - 1); 
            size_t currentCol = (nums[1] - 1); 

            if (currentRow < currentCol) {
                throw std::runtime_error("Matrix contains entries outside the lower triangular part and may not be symmetric");
            }

            IA[lineNum - 1] = currentRow;
            if (currentCol != nums[1]) {
                JA[currentCol] = lineNum - 1;
                currentCol++;
            }
            V[lineNum - 1] = nums[2];
        } else {
            size = nums[0];
            numEntries = nums[2];
            V.resize(numEntries);
            IA.resize(numEntries);
            JA.resize(size + 1);
        }
        nums.clear();
        lineNum++;
      }
    }
  }

  // Matrix-vector multiplication for symmetric matrices
  std::vector<double> multiplySymmetric(const std::vector<double> &vec) const {
    if (vec.size() != size) {
      throw std::runtime_error("Vector size does not match matrix dimensions.");
    }

    std::vector<double> result(size, 0.0);

    for (size_t j = 0; j < size; ++j) {
      for (size_t index = JA[j]; index < JA[j + 1]; ++index) {
        double val = V[index];
        size_t i = IA[index];

        result[i] += val * vec[j];

        if (i != j) {
          result[j] += val * vec[i];
        }
      }
    }

    return result;
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
  if (std::sscanf(argv[2], "%zu", &opts.iters) != 1 || opts.iters == 0)
    throw std::runtime_error("invalid parameter for iters");
  return opts;
}

} // namespace program_options

// Conjugate method
std::vector<double> conjugateGradient(const CCSmatrix& A, const std::vector<double>& b, size_t maxIters, double tol) {
    size_t n = b.size();
    std::vector<double> x(n, 0.0); 
    std::vector<double> r = b;  
    std::vector<double> p = r;  
    double r_norm, res;
    double r0_norm = sqrt(std::inner_product(b.begin(), b.end(), b.begin(), 0.0));
    double rTr_old = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
    std::vector<double> residuals;

    std::vector<double> x_star(n, 1.0);
    std::vector<double> errors;


    for (size_t k = 0; k < maxIters; ++k) {

        // 1. Step length alpha_k
        std::vector<double> Ap = A.multiplySymmetric(p);
        double alpha = rTr_old / std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0);

        // 2. Approximate solution x_k
        for (size_t i = 0; i < n; ++i)
            x[i] += alpha * p[i];

        // 3. Residual r_k
        for (size_t i = 0; i < n; ++i)
            r[i] -= alpha * Ap[i];

        // update residual vector
        // Calculate the residual norm
        double r_norm = sqrt(std::inner_product(r.begin(), r.end(), r.begin(), 0.0));
        res = r_norm / r0_norm;
        residuals.push_back (res);

        // 4. Beta_k  (improvement this step)
        double rTr_new = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
        double beta = rTr_new / rTr_old;

        // 5. Search direction p_k
        for (size_t i = 0; i < n; ++i)
            p[i] = r[i] + beta * p[i];

        rTr_old = rTr_new;

        // Compute error in A-norm
        std::vector<double> ek(n);
        for (size_t i = 0; i < n; ++i)
            ek[i] = x_star[i] - x[i]; // e_k = x_star - x

        std::vector<double> Aek = A.multiplySymmetric(ek);
        double err = sqrt(std::inner_product(ek.begin(), ek.end(), Aek.begin(), 0.0)); // sqrt(e_k^T * A * e_k)
        errors.push_back(err);

    }

    // Write residuals and errors to a text file
    std::ofstream outFile("residuals.txt");
    if (!outFile) {
        std::cerr << "Error: Unable to open residuals.txt for writing.\n";
    } else {
        outFile << "Iteration\tResidual\tError_A_norm\n";
        for (size_t i = 0; i < residuals.size(); ++i) {
            outFile << i + 1 << "\t" << residuals[i] << "\t" << errors[i] << "\n";
        }
        std::cout << "Results successfully written to residuals.txt\n";
    }
    outFile.close();
    
    double lastResidual = residuals.back();
    std::cout << "Computed residual-ratio: " << lastResidual << std::endl;

    double lastError= errors.back();
    std::cout << "Computed error in A-norm: " << lastError << std::endl;

    return x; 
}

int main(int argc, char *argv[]) try {
    // Parse program arguments
    auto opts = program_options::parse(argc, argv);
    opts.print();

    // Open and read the matrix file
    auto mat = CCSmatrix();
    std::ifstream matrixFile(opts.filename);
    if (!matrixFile.is_open()) {
        throw std::runtime_error("Unable to open the specified matrix file.");
    }
    mat.readMtxFile(matrixFile);

    // Test matrix-vector multiplication
    std::vector<double> x_star(mat.size, 1.0); // Example vector for multiplication

    std::vector<double> b = mat.multiplySymmetric(x_star);

    // Solve the system using Conjugate Gradient method
    std::vector<double> x = conjugateGradient(mat, b, opts.iters, 1e-6);

    return EXIT_SUCCESS;
} catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
}
