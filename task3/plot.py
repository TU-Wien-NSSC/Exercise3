import re
import matplotlib.pyplot as plt

# List of file paths
file_paths = ["SINX","COS2XINV","X4M5"]

# Colors or markers for each file
styles = ['o-', 's-', '^-']

# Plot Runtime vs No. of Threads for each file
plt.figure(figsize=(10,6))

for i, file_path in enumerate(file_paths):
    print(f"\nExtracting data from File: {file_path}")
    threads = []
    runtimes = []
    current_time = None

    # Extract data from file
    file = file_path + ".txt"
    print(f"\n{file}")
    with open(file, "r") as file:
        for line in file:
            # Match "Time: <value>"
            time_match = re.search(r"Runtime:\s*([\d.]+)", line)
            if time_match:
                current_time = float(time_match.group(1))
            
            # Match "Threads: <value>"
            thread_match = re.search(r"Threads:\s(\d+)", line)
            if thread_match and current_time is not None:
                threads.append(int(thread_match.group(1)))
                runtimes.append(current_time)
                current_time = None

    
    for rt, th in zip(runtimes, threads):
        print(f"Threads: {th}, Runtime: {rt}")
    plt.plot(threads, runtimes, styles[i%len(styles)],label=f"{file_path}")

plt.xlabel("Number of Threads")
plt.ylabel("Runtime (seconds)")
plt.title("Runtime vs No. of Threads for Monte Carlo Integration")
plt.legend()
plt.grid(True)

# Save the plot as a PNG file
output_file = "runtime_vs_threads.png"
plt.savefig(output_file)

plt.show()