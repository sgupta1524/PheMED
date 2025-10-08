import subprocess
import numpy as np
import pandas as pd
import re
import json
import argparse
parser = argparse.ArgumentParser()

def run_phemed_with_conditions(initial_conditions, optimisation_method, output_dir, sum_stats, n_studies):
    if isinstance(initial_conditions, (list, tuple)):
        init_str = ",".join(str(v) for v in initial_conditions)
    else:
        raise ValueError("initial_conditions must be a list or tuple of numbers")
    command = [
        "python", "phemed.py",
        "--sum_stats", str(sum_stats),
        "--n_studies", str(n_studies),
        "--optimisation_init", init_str,
        "--optimizer_method", str(optimisation_method),
        "--out", str(output_dir),
        "--compute_cis", str(False)
    ]
    print("Running:", " ".join(command))
    
    process = subprocess.Popen(command, stdout=subprocess.PIPE, 
                              stderr=subprocess.STDOUT, text=True)
    
    dilution_vals = None
    for line in iter(process.stdout.readline, ''):
        if "Effective dilution values are :" in line:
            match = re.search(r"\[(.*?)\]", line)
            if match:
                dilution_vals = [float(x.strip()) for x in match.group(1).split(',')]
                break
    
    process.wait()
    return dilution_vals

# write s function to check if diltuion values above are withon +- 0.5 of each other
def check_dilution_stability(dilution_vals1, dilution_vals2, tolerance):
    print(dilution_vals1, dilution_vals2)
    if dilution_vals1 is None or dilution_vals2 is None:
        return "invalid inputs"
    if all(abs(a - b) <= tolerance for a, b in zip(dilution_vals1, dilution_vals2)):
        return "within tolerance limit"
    else:
        return "not all within tolerance limit"

def main():
    parser = argparse.ArgumentParser(description="Run PheMED stability analysis.")
    parser.add_argument("--sum_stats", required=True, help="Path to summary statistics file.")
    parser.add_argument("--n_studies", required=True, type=int, help="Number of studies.")
    parser.add_argument("--output_dir", required=True, help="Directory to store output files.")
    parser.add_argument("--initial_conditions_file", required=True, default="initial_conditions_list.json", 
                        help="JSON file with list of initial conditions. Please see the example file " \
                        "initial_conditions_list.json")
    parser.add_argument("--tolerance", type=float, default=0.5, help="Tolerance for stability check.")
    args = parser.parse_args()

    sum_stats = args.sum_stats
    n_studies = args.n_studies
    output_dir = args.output_dir
    init_file = args.initial_conditions_file
    tolerance = args.tolerance
    # run run_phemed_with_conditions with five different initial conditions and optimisation methods
    with open(init_file, "r") as file:
        initial_conditions_list = json.load(file)
        optimisation_methods = [
            "Nelder-Mead",
            "BFGS",
            "L-BFGS-B",
            "SLSQP"
        ]
        dilution_results = []
        for initial_conditions in initial_conditions_list:
            for optimisation_method in optimisation_methods:
                dilution_vals = run_phemed_with_conditions(initial_conditions, optimisation_method, output_dir, sum_stats, n_studies)
                dilution_results.append((initial_conditions, optimisation_method, dilution_vals))
                #run check_dilution_stability on results pairwise
        for i in range(len(dilution_results)):
            for j in range(i + 1, len(dilution_results)):
                cond1, method1, vals1 = dilution_results[i]
                cond2, method2, vals2 = dilution_results[j]
                stable = check_dilution_stability(vals1, vals2, tolerance)
                print(f"Stability between {method1} ({cond1}) and {method2} ({cond2}): {stable}")
                if stable == "not all within tolerance limit":
                    print("Significant instability detected in phemed, stopping further checks.")
                    break
                else:
                    print("No significant instability detected in phemed.")
                    continue

if __name__ == "__main__":
    main()