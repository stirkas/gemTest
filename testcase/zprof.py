import sys
import matplotlib.pyplot as plt
import numpy as np

def read_input(file_path):
    var1_values = []
    var2_values = []
    var3_values = []

    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 3):
            # Assuming timestamp is at the beginning of each line
            timestamp, *var1 = map(float, lines[i].split())
            var1_values.append((timestamp, var1))

            _, *var2 = map(float, lines[i + 1].split())
            var2_values.append((timestamp, var2))

            _, *var3 = map(float, lines[i+1].split())
            var3_values.append((timestamp,var3))

    return var1_values, var2_values, var3_values

def plot_at_timestamp(var1_values, var2_values, var3_values, timestamp_index):
    timestamps_var1, values_var1 = zip(*var1_values)
    timestamps_var2, values_var2 = zip(*var2_values)
    timestamps_var3, values_var3 = zip(*var3_values)

    timestamp_at_index = timestamps_var1[timestamp_index]
    xarr = np.linspace(0.2,0.8,128)
    intFac = float(int(timestamp_index) - 4)
    vals1 = values_var1[timestamp_index]
    vals2 = values_var2[timestamp_index]
    vals3 = values_var3[timestamp_index]
    vals1 = [x/intFac for x in vals1]
    vals2 = [x/intFac for x in vals2]
    vals3 = [x/intFac for x in vals3]
    plt.plot(xarr, vals1, label='$\\delta n_{e,SM}$')
    plt.plot(xarr, vals2, label='$\\delta T_{e,SM}$')
    plt.plot(xarr, vals3, label='$\\delta E_{e,SM}$')
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

    var1_values, var2_values, var3_values  = read_input(file_path)
    
    if 0 <= timestamp_index < len(var1_values):
        plot_at_timestamp(var1_values, var2_values, var3_values, timestamp_index)
    else:
        print(f"Invalid timestamp index. It should be between 0 and {len(var1_values) - 1}.")

