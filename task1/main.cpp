#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <numeric>


struct CCSMatrix {
    std::vector<double> V;
    std::vector<int> IA;
    std::vector<int> JA;
    size_t size;

    void readMtxFile(const std::string &filename) {
        std::ifstream file(filename);

        if (!file.is_open()) {
            throw std::runtime_error("Error: Unable to open file.");
        }

        std::string line;
        std::getline(file, line); 
        int M, N, L; 
        file >> N >> M >> L;
        size = N;

        // Initialize column pointers
        JA.resize(N + 1, 0);

        int currentRow, currentCol;
        double currentVal;
        int lastIndex = -1;

        // Read the entries
        while (L--) {
            file >> currentRow >> currentCol >> currentVal;
            currentRow--; // Adjusting to 0-based indexing
            currentCol--;

            // Fill the values and the rows vectors
            V.push_back(currentVal);
            IA.push_back(currentRow);

            // Fill the columns vector
            if (currentCol != lastIndex) {
                JA[currentCol] = V.size() - 1;
                lastIndex = currentCol;
            }
        }

        // Finalize column pointers
        JA.back() = V.size();
    }
};

// Matrix-vector multiplication for symmetric matrices
std::vector<double> multiplySymmetric(const CCSMatrix& matrix, const std::vector<double>& vec) {
    if (vec.size() != matrix.size) {
        throw std::runtime_error("Vector size does not match matrix dimensions.");
    }

    std::vector<double> result(matrix.size, 0.0); // Initialize result vector with zeros

    for (size_t j = 0; j < matrix.size; ++j) {
        for (size_t k = matrix.JA[j]; k < matrix.JA[j + 1]; ++k) {
            size_t i = matrix.IA[k]; 

            result[i] += matrix.V[k] * vec[j];

            if (i != j) {
                result[j] += matrix.V[k] * vec[i];
            }
        }
    }

    return result;
}

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

} 

// Conjugate method
std::vector<double> conjugateGradient(const CCSMatrix& A, const std::vector<double>& b, const std::vector<double>& x_star, size_t maxIters) {
    size_t n = b.size();
    std::vector<double> x(n, 0.0); 
    std::vector<double> r = b;  
    std::vector<double> p = r;  
    double r_norm, res, err;
    double r0_norm = sqrt(std::inner_product(b.begin(), b.end(), b.begin(), 0.0));
    double rTr_old = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
    



    // Write residuals and errors to a text file
    std::ofstream outFile("residuals_and_errors.txt");
    if (!outFile) {
        std::cerr << "Error: Unable to open residuals.txt for writing.\n";
    } else {

            outFile << "Iteration\tResidual\tError_A_norm\n";
            
            for (size_t k = 0; k < maxIters; ++k) {

            // 1. Step length alpha_k
            std::vector<double> Ap = multiplySymmetric(A, p);
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

            std::vector<double> Aek = multiplySymmetric(A, ek);
            err = sqrt(std::inner_product(ek.begin(), ek.end(), Aek.begin(), 0.0)); // sqrt(e_k^T * A * e_k)

            outFile << k + 1 << "\t" << res << "\t" << err << "\n";

        }
        double lastResidual = res;
        std::cout << "Last computed residual-ratio: " << lastResidual << std::endl;

        double lastError= err;
        std::cout << "Last computed error in A-norm: " << lastError << std::endl;

    }
    outFile.close();
    
    

    return x; 
}

// Main function to test the CCS function (Takes as input the name of the .mtx file and the number of iteration for the CG method)
int main(int argc, char* argv[]) {
    // Parse program arguments
    auto opts = program_options::parse(argc, argv);
    opts.print();

    // Create the matrix and read from file
    CCSMatrix mat;
    mat.readMtxFile(opts.filename);

    // Test matrix-vector multiplication
    std::vector<double> x_star(mat.size, 1.0);    
    std::vector<double> b = multiplySymmetric(mat, x_star);

    // Solve the system using Conjugate Gradient method
    std::vector<double> x = conjugateGradient(mat, b, x_star, opts.iters);

    return 0;
}
