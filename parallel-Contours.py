import subprocess
from concurrent.futures import ThreadPoolExecutor
import sys
import os

def run_script(args):
    subprocess.run(["python", "-m", "src.data_extraction.DataCollector_Deterministic", str(args[0]), str(args[1])])


def main(path, variable1, variable2):
    arguments = [(os.path.join(path, "inputs", f'{variable1}_{variable2}.{a}_{b}.json'), path) for a, b in split_digits(20)]

    with ThreadPoolExecutor(max_workers=100) as executor:
        executor.map(run_script, arguments)

def split_digits(stop):
    digits = []
    for i in range(stop):
        b = i % 10
        a = i // 10
        digits.append((a,b))
    return digits

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'Usage: python3 {sys.argv[0]} <folder_path> <mode>')
        sys.exit(1)

    input_names = ["delay", "relation", "R0"]

    folder_path = sys.argv[1]
    mode = int(sys.argv[2])
    var1 = input_names[mode % 3]
    var2 = input_names[(mode + 1) % 3]
    main(folder_path, var1, var2)

