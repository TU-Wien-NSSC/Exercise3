import os
import re

# Directory containing the output files
output_dir = "results"

# List of threads and scheduling policies
threads = [1, 2, 4, 8, 10, 20, 40]
policies = ["dynamic", "static", "static,1"]

# Regular expression to extract runtime
runtime_pattern = re.compile(r"Total runtime: ([\d.]+) seconds")

# Dictionary to store runtimes
runtimes = {policy: {} for policy in policies}

# Debugging: Ensure the output directory exists
if not os.path.exists(output_dir):
    print(f"Error: Directory '{output_dir}' does not exist.")
else:
    # Loop through all files to extract runtime
    for t in threads:
        for policy in policies:
            filename = f"solver_output_{t}_{policy}.txt"
            filepath = os.path.join(output_dir, filename)
            if os.path.exists(filepath):
                with open(filepath, 'r') as file:
                    content = file.read()
                    match = runtime_pattern.search(content)
                    if match:
                        runtime = float(match.group(1))
                        runtimes[policy][t] = runtime
            else:
                print(f"Warning: File '{filename}' does not exist in '{output_dir}'.")

# Save data to a file
output_file = "runtimes.csv"
with open(output_file, 'w') as f:
    f.write("Threads,Policy,Runtime\n")
    for policy in policies:
        for t in threads:
            if t in runtimes[policy]:
                f.write(f"{t},{policy},{runtimes[policy][t]}\n")
print(f"Runtimes saved to {output_file}")
