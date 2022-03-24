import numpy as np
import matplotlib.pyplot as plt
import sys


path_to_dir = sys.argv[1]
PLOT = True


def read_times(path_to_file):
    with open(path_to_file, "r") as myfile:
        times = myfile.read()[:-1].split(",")
    return list(float(n) for n in times)


def plot(lTrue, lFalse, title, avgTrue, avgFalse):
    fig = plt.figure(figsize=(9,6))
    plt.plot(lTrue[0], "-x", label="Ids in C, logged, avg %fs"%(avgTrue[0]))
    plt.plot(lFalse[0], "-o", label="Ids in C, not logged, avg %fs"%(avgFalse[0]))
    plt.plot(lTrue[1], "-*", label="Ids in Python, logged, avg %fs"%(avgTrue[1]))
    plt.plot(lFalse[1], ".", label="Ids in Python, not logged, avg %fs"%(avgFalse[1]))
    plt.xlabel("Number of run")
    plt.ylabel("Time [s]")
    plt.legend()
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path_to_dir + "results.png")
    plt.show()

# collect and compute data
avg_diff_compile = np.zeros(2)
avg_diff_compute = np.zeros(2)
perc_diff_compute = np.zeros(2)
avg_logTrue_compute = np.zeros(2)
avg_logFalse_compute = np.zeros(2)
logTrue_compute = []
logFalse_compute = []
for i, commit in enumerate(["old", "new"]):
    run = "compile" 
    path_to_file = f"{path_to_dir}time_log%s_{run}_{commit}.txt"%True
    times_logTrue = read_times(path_to_file)
    path_to_file = f"{path_to_dir}time_log%s_{run}_{commit}.txt"%False
    times_logFalse = read_times(path_to_file)

    avg_diff_compile[i] = np.average([float(n)-float(m) for n,m in zip(times_logTrue, times_logFalse)])
    
    run = "compute"
    path_to_file = f"{path_to_dir}time_log%s_{run}_{commit}.txt"%True
    times_logTrue = read_times(path_to_file)
    path_to_file = f"{path_to_dir}time_log%s_{run}_{commit}.txt"%False
    times_logFalse = read_times(path_to_file)

    avg_logTrue_compute[i] = np.average([n for n in times_logTrue])
    avg_logFalse_compute[i] = np.average([n for n in times_logFalse])
    avg_diff_compute[i] = avg_logTrue_compute[i] - np.average([n for n in times_logFalse])
    perc_diff_compute[i] = 100*(1 - np.average([n for n in times_logFalse])/np.average([n for n in times_logTrue]))
    logTrue_compute += [times_logTrue]
    logFalse_compute += [times_logFalse]

perc_speedup_compute = 100*(1 - avg_logTrue_compute[1]/avg_logTrue_compute[0])
abs_speedup_compute = avg_logTrue_compute[0]-avg_logTrue_compute[1]

# plot and print results
if PLOT:
    title = (f"""How much faster is run with no logging compared to with logging when ids set in C? \nabsolut {avg_diff_compute[0]:3.4}, percent  {perc_diff_compute[0]:3.4}\n"""
            f"""How much faster is run with no logging compared to with logging when ids set in python?\nabsolut {avg_diff_compute[1]:3.4}, percent {perc_diff_compute[1]:3.4}\n""" 
            f"""How much faster is run with logging when ids set in python compared to when it is set in C?\nabsolut {abs_speedup_compute:3.4}, percent {perc_speedup_compute:3.4}\n""")
    plot(logTrue_compute, logFalse_compute, title, avg_logTrue_compute, avg_logFalse_compute)

print("\n -------RESULTS---------\n")
print("Compile abs time averaged 0=old 1=new", avg_diff_compile)
print("Compute abs time averaged 0=old 1=new", avg_diff_compute)
print("How much faster is no logging in each commit (0=old 1=new) in percent?", perc_diff_compute)
print("How much faster is logging between the commits in percent?", perc_speedup_compute)
print("How much faster is logging between the commits absolute?", abs_speedup_compute)