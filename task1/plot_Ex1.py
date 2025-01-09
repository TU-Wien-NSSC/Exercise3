import matplotlib.pyplot as plt

# Read data from the residuals.txt file
iterations = []
residuals = []
errors = []

# Define down-sampling step
step = 100

with open("residuals_and_errors.txt", "r") as file:
    # Skip the header
    next(file)
    for i, line in enumerate(file):
        if i % step == 0:  
            parts = line.strip().split("\t")
            iterations.append(int(parts[0]))
            residuals.append(float(parts[1]))
            errors.append(float(parts[2]))

# Plot residuals
plt.figure(figsize=(10, 5))
plt.semilogy(iterations, residuals, label="Residuals", color="blue")
plt.xlabel("Iteration")
plt.ylabel("Residual")
plt.title("Convergence of Residuals")
plt.grid(True)
plt.savefig("PlotResidual.png")

# Plot errors in A-norm
plt.figure(figsize=(10, 5))
plt.semilogy(iterations, errors, label="Error (A-norm)", color="red")
plt.xlabel("Iteration")
plt.ylabel("Error (A-norm)")
plt.title("Convergence of Errors in A-Norm")
plt.grid(True)
plt.savefig("PlotError.png")
plt.show()
