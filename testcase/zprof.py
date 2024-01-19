import sys
import matplotlib.pyplot as plt
import numpy as np

def read_input(file_path):
    var1_values = []
    var2_values = []

    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 2):
            # Assuming timestamp is at the beginning of each line
            timestamp, *var1 = map(float, lines[i].split())
            var1_values.append((timestamp, var1))

            _, *var2 = map(float, lines[i + 1].split())
            var2_values.append((timestamp, var2))

    return var1_values, var2_values

def plot_at_timestamp(var1_values, var2_values, timestamp_index):
    timestamps_var1, values_var1 = zip(*var1_values)
    timestamps_var2, values_var2 = zip(*var2_values)

    timestamp_at_index = timestamps_var1[timestamp_index]
    xarr = 0.2 + np.arange(1,129)*(0.6/127)
    plt.plot(xarr, values_var1[timestamp_index], label='$\\delta n_{e,SM}$')
    plt.plot(xarr, values_var2[timestamp_index], label='$\\delta T_{e,SM}$')
    plt.xlabel('r/a')
    plt.legend()
    plt.title(f't = {timestamp_at_index}')
    plt.grid()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <file_path> <timestamp_index>")
        sys.exit(1)

    file_path = sys.argv[1]
    timestamp_index = int(sys.argv[2])

    var1_values, var2_values = read_input(file_path)
    
    if 0 <= timestamp_index < len(var1_values):
        plot_at_timestamp(var1_values, var2_values, timestamp_index)
    else:
        print(f"Invalid timestamp index. It should be between 0 and {len(var1_values) - 1}.")

