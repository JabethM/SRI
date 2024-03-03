import subprocess
from concurrent.futures import ThreadPoolExecutor
import sys
import os


def run_script(args):
    subprocess.run(
        ["python", "-m", "src.data_extraction.DataCollector_Deterministic", str(args[0]), str(args[1]), str(args[2])])


def main(path, variable1, variable2, batch, name_suffix):
    arguments = [(os.path.join(path, "inputs", f'{variable1}_{variable2}-{batch}{b}.json'),
                  path, f"{name_suffix}.{batch}.{b}") for b in range(10)]

    with ThreadPoolExecutor(max_workers=10) as executor:
        executor.map(run_script, arguments)


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print(f'Usage: python3 {sys.argv[0]} <folder_path> <mode> <batch> <suffix>')
        sys.exit(1)

    input_names = ["delay", "relation", "R0"]

    folder_path = sys.argv[1]
    mode = int(sys.argv[2])
    var1 = input_names[mode % 3]
    var2 = input_names[(mode + 1) % 3]
    data_batch = int(sys.argv[3])
    suffix = sys.argv[4]
    main(folder_path, var1, var2, data_batch, suffix)
