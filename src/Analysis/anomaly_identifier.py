import os
import csv


def find_matching_files(directory):
    anomalous = []
    single_var = []
    for filename in os.listdir(directory):
        if filename.startswith("dataFile_") and filename.endswith(".csv"):
            filepath = os.path.join(directory, filename)
            if is_matching_file(filepath) == 1:
                anomalous.append(filename)
            elif is_matching_file(filepath) == 2:
                single_var.append(filename)
            else:
                continue
    return anomalous, single_var


def is_matching_file(filepath):
    with open(filepath, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            if row[0] == "Peak_Infections":
                times = list(map(float, row[1].split(',')))
                if times[0] <= 5 and times[1] <= 5:
                    return 1 # Anomaly
                elif times[0] > 5 >= times[1]:
                    return 2  # Single Var run
                else:
                    return 0
    return False


directory = r"D:\MPhys\Github\First-Take\SRI-data\Watts\Rel.vs.R0\LowDelay"
a, s = find_matching_files(directory)
print("- Single Files: ")
for file in s:
    print(file)
print("- Anomalous files:")
for file in a:
    print(file)
