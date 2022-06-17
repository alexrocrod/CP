import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

save = True

dir_name = os.getcwd() + "\cp_harrisDetector"

df = pd.read_excel("host_times.xlsx")

plt.close("all")

image_names = ["chess", "chessBig", "chessL", "chessRotate1", "house"]
summary_results_df = pd.DataFrame([], index=image_names, columns=["Host", "OpenMP", "CUDA Global Memory", "CUDA Constant Memory"])

for file in os.listdir(dir_name):
    if file.endswith(".txt") and "reference" not in file:

        split_file_name = file.split("_")

        if "CUDA" in file:
            module_name = split_file_name[0]
            memory_type = split_file_name[1]
            image_name = split_file_name[4].split(".")[0]
        else:
            module_name = split_file_name[0]
            title = module_name
            image_name = split_file_name[3].split(".")[0]

        host_time = df[df["Image"] == image_name][module_name].values[0]

        with open(dir_name + f"\{file}") as f:
            lines = f.readlines()

        times = [float(l.strip("\n")) for l in lines]
        mean_time = np.mean(times)

        host_time = round(host_time, 3)
        mean_time = round(mean_time, 3)

        if module_name == "OpenMP":
            summary_results_df.loc[image_name, "OpenMP"] = mean_time
            summary_results_df.loc[image_name, "Host"] = host_time
        else:
            summary_results_df.loc[image_name, f"CUDA {memory_type} Memory"] = mean_time
            summary_results_df.loc[image_name, "Host"] = host_time
        
        plt.figure()
        plt.plot(np.arange(1, 101), times)
        plt.axhline(mean_time, c="r", ls = "--",label=f"Tempo médio: {mean_time:.2f} (ms)")
        plt.xlabel("run")
        plt.ylabel("tempo")
        if "CUDA" in file:
            plt.title(f"Tempos em função das runs da imagem {image_name}.pgm com {module_name} {memory_type} memory\nTempo host: {host_time:.2f} (ms)")
        else:
            plt.title(f"Tempos em função das runs da imagem {image_name}.pgm com {module_name}\nTempo host: {host_time:.2f} (ms)")
        plt.legend()

        if save:
            if "CUDA" in file:
                plt.savefig(f"figures/{module_name}_{memory_type}_{image_name}_plot_times.jpeg", bbox_inches="tight")
            else:
                plt.savefig(f"figures/{module_name}_{image_name}_plot_times.jpeg", bbox_inches="tight")

#plt.show()

summary_results_df.to_csv("summary_of_time_results.csv")

speed_up_df = pd.DataFrame([], index=summary_results_df.index, columns=["OpenMP", "CUDA Global Memory", "CUDA Constant Memory"])


speed_up_df["OpenMP"] = summary_results_df["Host"] / summary_results_df["OpenMP"]
speed_up_df["CUDA Global Memory"] = summary_results_df["Host"] / summary_results_df["CUDA Global Memory"]
speed_up_df["CUDA Constant Memory"] = summary_results_df["Host"] / summary_results_df["CUDA Constant Memory"]

speed_up_df.to_csv("summary_speed_up_results.csv")

res_df = pd.read_excel("resOpenMP.xlsx")
#print(res_df)

n_threads = np.arange(1, 25)
module = res_df["module"].values[0]

for img in res_df["image"].unique():
    times = res_df[res_df["image"] == img].time.values
    optimal_n_thread_ind = np.argmin(times)
    plt.figure()
    plt.plot(n_threads, times)
    plt.plot(n_threads[optimal_n_thread_ind], times[optimal_n_thread_ind], "r.", label = f"Número ótimo de threads: {n_threads[optimal_n_thread_ind]}")
    plt.xlabel("nº de threads")
    plt.ylabel("tempo médio (ms)")
    plt.title(f"Tempo médio em função do nº de threads\nda imagem {img}.pgm com {module}")
    plt.xticks(n_threads)
    plt.legend()
    if save:
        plt.savefig(f"figures/{module_name}_{img}_plot_n_threads.jpeg", bbox_inches="tight")

#plt.show()