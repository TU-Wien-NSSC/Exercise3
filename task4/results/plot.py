import csv
import matplotlib.pyplot as plt

file_name = "main"
#file_name = "main_ex2"

file_path = f"results_{file_name}.txt"  

# Initialize a dictionary to store data by policy
data = {"dynamic": [], "static": [], "static_1": []}

# Read the CSV file
with open(file_path, newline='') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header row
    for row in reader:
        threads, policy, runtime = int(row[0]), row[1], float(row[2])
        data[policy].append((threads, runtime))

# Calculate speedup for each policy
speedup_data = {}
for policy in data:
    # Get the runtime for 1 thread (serial execution)
    serial_runtime = next(runtime for threads, runtime in data[policy] if threads == 1)
    
    # Calculate speedup for each entry
    speedup_data[policy] = [
        (threads, serial_runtime / runtime) for threads, runtime in data[policy]
    ]

# Plot the speedup for each policy
plt.figure(figsize=(8, 6))
for policy, speedups in speedup_data.items():
    # Separate the thread counts and speedup values
    threads = [entry[0] for entry in speedups]
    speedup = [entry[1] for entry in speedups]
    plt.plot(threads, speedup, marker="o", label=policy.replace("_", ","))

# Add labels, legend, and title
plt.xlabel("Number of Threads")
plt.ylabel("Speedup")
plt.title("Speedup vs Threads for Different Scheduling Policies")
plt.legend(title="Scheduling Policy")
plt.grid(True)

# Save and show the plot
plt.savefig(f"speedup_plot_{file_name}.png")
plt.show()
