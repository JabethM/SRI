import multiprocessing

import sys
import os


def run_script(args):
    os.system(f"python3 {os.path.join('src', 'parallel-Contours-Batch.py')}"
              f" {str(args[0])} {str(args[1])} {str(args[2])} {str(args[3])}")


def main(path, mode, batch, repeats):
    num_of_repeats = repeats
    arguments = [(path, mode, batch, b) for b in range(num_of_repeats)]

    pool = multiprocessing.Pool(processes=num_of_repeats)
    pool.map(run_script, arguments)

    pool.close()
    pool.join()


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print(f'Usage: python3 {sys.argv[0]} <folder_path> <mode> <batch> <repeats>')
        sys.exit(1)

    folder_path = sys.argv[1]
    contour_mode = int(sys.argv[2])
    data_batch = int(sys.argv[3])
    runs = int(sys.argv[4])
    main(folder_path, contour_mode, data_batch, runs)
