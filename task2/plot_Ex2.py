import matplotlib.pyplot as plt
import numpy as np

# Define a function to read the data from a file
def read_residuals(filename):
    iterations = []
    residuals = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith("Iteration"):  # Skip header or empty lines
                iteration, residual = line.split()
                iterations.append(int(iteration))
                residuals.append(float(residual))
    return iterations, residuals

# Read data from both files
iterations_noprecond, residuals_noprecond = read_residuals("residual_noprecond.txt")
iterations_withprecond, residuals_withprecond = read_residuals("residual_withprecond.txt")

# Plot the data
plt.figure(figsize=(10, 6))
plt.semilogy(iterations_noprecond, residuals_noprecond, label="No Precondition", color='blue')
plt.semilogy(iterations_withprecond, residuals_withprecond, label="With Diagonal Precond.", color='green')

# Customize the plot
plt.title("Convergence of Residuals")
plt.xlabel("Iteration")
plt.ylabel("Residual")
plt.legend()
plt.savefig("PlotResPrecond.png")

# Show the plot
plt.show()
