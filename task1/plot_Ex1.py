import matplotlib.pyplot as plt

# Read data from the residuals.txt file
iterations = []
residuals = []
errors = []

with open("residuals.txt", "r") as file:
    # Skip the header
    next(file)
    for line in file:
        parts = line.strip().split("\t")
        iterations.append(int(parts[0]))
        residuals.append(float(parts[1]))
        errors.append(float(parts[2]))

# Plot residuals
plt.figure(figsize=(10, 6))
plt.plot(iterations, residuals, label="Residuals", marker="o", linestyle="-", color="blue")
plt.xlabel("Iteration")
plt.ylabel("Residual")
plt.title("Convergence of Residuals")
plt.grid(True)
plt.savefig("residuals_plot.png")
plt.show()

# Plot errors in A-norm
plt.figure(figsize=(10, 6))
plt.plot(iterations, errors, label="Error (A-norm)", marker="x", linestyle="--", color="red")
plt.xlabel("Iteration")
plt.ylabel("Error (A-norm)")
plt.title("Convergence of Errors in A-Norm")
plt.grid(True)
plt.savefig("errors_plot.png")
plt.show()
