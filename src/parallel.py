import subprocess
from concurrent.futures import ThreadPoolExecutor
import sys
import os


def run_script(args):
    subprocess.run(["python", "-m", "src.data_extraction.DataCollector_Deterministic", str(args[0]), str(args[1])])


def main(path):
    arguments = [(os.path.join(path, f'transmission_rate_{i}.json'), os.path.join(path, "data")) for i in range(20)]

    with ThreadPoolExecutor(max_workers=20) as executor:
        executor.map(run_script, arguments)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f'Usage: python3 {sys.argv[0]} <folder_path>')
        sys.exit(1)

    folder_path = sys.argv[1]
    main(folder_path)

