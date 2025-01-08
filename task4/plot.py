import matplotlib.pyplot as plt
import csv
from collections import defaultdict

# Input CSV file
input_file = "runtimes.csv"

# List of threads
threads = [1, 2, 4, 8, 10, 20, 40]

# Dictionary to store runtimes
runtimes = defaultdict(dict)

# Read the runtimes from the CSV file
with open(input_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        thread_count = int(row["Threads"])
        policy = row["Policy"]
        runtime = float(row["Runtime"])
        runtimes[policy][thread_count] = runtime

# Plot the speedup vs threads
plt.figure(figsize=(10, 6))
for policy, policy_runtimes in runtimes.items():
    if 1 not in policy_runtimes:
        print(f"Skipping policy '{policy}' as it does not have a baseline runtime for 1 thread.")
        continue
    baseline_runtime = policy_runtimes[1]
    speedup = []
    valid_threads = []
    for t in threads:
        if t in policy_runtimes:
            valid_threads.append(t)
            speedup.append(baseline_runtime / policy_runtimes[t])
    if speedup:
        plt.plot(valid_threads, speedup, marker='o', label=f"Policy: {policy}")

# Formatting the plot
plt.title("Speedup vs Number of Threads")
plt.xlabel("Number of Threads")
plt.ylabel("Speedup")
plt.xscale("log", base=2)
plt.xticks(threads)
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend()
plt.tight_layout()

# Save the plot to a file
plot_filename = "speedup_plot.png"
plt.savefig(plot_filename)
print(f"Speedup plot saved to {plot_filename}")

# Show the plot
plt.show()
