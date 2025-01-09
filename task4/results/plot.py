import csv
import matplotlib.pyplot as plt

#change following to create different plot

file_name = "main"
#file_name = "main_ex2"

file_path = f"results_{file_name}.txt"  
data = {"dynamic": [], "static": [], "static_1": []}


with open(file_path, newline='') as file:
    reader = csv.reader(file)
    next(reader)  #skip first row
    for row in reader:
        threads, policy, runtime = int(row[0]), row[1], float(row[2])
        data[policy].append((threads, runtime))


speedup_data = {}
efficiency_data = {}
for policy in data:
    #get the runtime for 1 thread (serial execution)
    serial_runtime = next(runtime for threads, runtime in data[policy] if threads == 1)
    
    #calculate speedup for each entry relative to the runtime with one thread
    speedup_data[policy] = [
        (threads, serial_runtime / runtime) for threads, runtime in data[policy]
    ]

    efficiency_data[policy] = [
        (threads, (serial_runtime / runtime) / threads) for threads, runtime in data[policy]
    ]


#Speedup Plot
plt.figure(figsize=(10, 7))
for policy, speedups in speedup_data.items():
    threads = [entry[0] for entry in speedups]
    speedup = [entry[1] for entry in speedups]
    plt.plot(threads, speedup, marker="o", label=policy.replace("_", ","))

plt.xlabel("Number of Threads")
plt.ylabel("Speedup")
plt.title("Scaling Speedup")
plt.legend(title="Scheduling Policy")
plt.grid(True)

plt.savefig(f"speedup_plot_{file_name}.png")
plt.show()


#Runtime Plot
plt.figure(figsize=(10, 7))
for policy, runtimes in data.items():
    threads = [entry[0] for entry in runtimes]
    runtime = [entry[1] for entry in runtimes]
    plt.plot(threads, runtime, marker="o", label=policy.replace("_", ","))

plt.xlabel("Number of Threads")
plt.ylabel("Runtime (seconds)")
plt.title("Runtime vs Threads")
plt.legend(title="Scheduling Policy")
plt.grid(True)

plt.savefig(f"runtime_plot_{file_name}.png")
plt.show()

#efficiency Plot
plt.figure(figsize=(10, 7))
for policy, efficiencies in efficiency_data.items():
    threads = [entry[0] for entry in efficiencies]
    efficiency = [entry[1] for entry in efficiencies]
    plt.plot(threads, efficiency, marker="o", label=policy.replace("_", ","))

plt.xlabel("Number of Threads")
plt.ylabel("Efficiency")
plt.title("Scaling Efficiency")
plt.legend(title="Scheduling Policy")
plt.grid(True)

plt.savefig(f"efficiency_plot_{file_name}.png")
plt.show()
