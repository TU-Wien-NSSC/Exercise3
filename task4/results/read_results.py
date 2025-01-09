import os
import re

policies = ["dynamic", "static", "static,1"]
main_files = {policy: f"task4_results_{policy}_main.txt" for policy in policies}
main_ex2_files = {policy: f"task4_results_{policy}_main_ex2.txt" for policy in policies}

output_main = "results_main.txt"
output_main_ex2 = "results_main_ex2.txt"

#get runtime from file
def extract_runtime(file_path, pattern):
    data = []
    if os.path.exists(file_path):
        with open(file_path, "r") as file:
            content = file.read()
            matches = pattern.findall(content)
            for threads, runtime in matches:
                data.append((int(threads), float(runtime)))
    else:
        print(f"File not found: {file_path}")
    return data

#pattern for runtime extraction
main_pattern = re.compile(r"Threads:\s+(\d+).*?TotalRuntime:\s+([\d.]+)\s+seconds", re.DOTALL)
main_ex2_pattern = re.compile(r"Threads:\s+(\d+).*?Runtime:\s+([\d.]+)", re.DOTALL)

#save results for main files
with open(output_main, "w") as f_main:
    f_main.write("Threads,Policy,Runtime\n")
    for policy, file_path in main_files.items():
        data = extract_runtime(file_path, main_pattern)
        for threads, runtime in data:
            if policy == "static,1":
                policy = "static_1"
            f_main.write(f"{threads},{policy},{runtime}\n")

# save results for main_ex2 files
with open(output_main_ex2, "w") as f_main_ex2:
    f_main_ex2.write("Threads,Policy,Runtime\n")
    for policy, file_path in main_ex2_files.items():
        data = extract_runtime(file_path, main_ex2_pattern)
        for threads, runtime in data:
            if policy == "static,1":
                policy = "static_1"
            f_main_ex2.write(f"{threads},{policy},{runtime}\n")

print(f"Results saved to {output_main} and {output_main_ex2}")
