import os
import sys


def parse_filename(filename):
    """Parse numbers from filename."""
    parts = filename.split('_')[-1].split('.')[:3]
    return tuple(map(int, parts))


def get_anomalies(b, c, filename):
    """Extract existing numbers from the file after '- Anomalous files:'."""
    existing_numbers = []
    start_reading = False # Set True to exclude single files
    with open(filename, 'r') as file:
        for line in file:
            if line.strip() == '- Single Files:':
                continue
            if line.strip() == '- Anomalous files:':
                start_reading = True # Set False to include Anomalous files
                continue
            if start_reading and line.strip():
                number = parse_filename(line.strip())
                if number[1] == b and number[2] == c:
                    existing_numbers.append(number[0])
    return existing_numbers


def filter_numbers(b, c, filename):
    """Filter out numbers based on file content."""
    existing_numbers = get_anomalies(int(b), int(c), filename)
    filtered_numbers = [num for num in range(10) if num not in existing_numbers]
    return filtered_numbers


def generate_file_paths(input_file, anom_file):
    file_paths = []
    directory, filename = os.path.split(input_file)
    base_name, extension = os.path.splitext(filename)
    version_parts = base_name.split('_')[1].split('.')
    repeats = filter_numbers(version_parts[1], version_parts[2], anom_file)
    #repeats = range(10)
    for i in repeats:
        version_parts[0] = str(i)
        new_filename = base_name.split('_')[0] + '_' + '.'.join(version_parts) + extension
        file_paths.append(os.path.join(directory, new_filename))

    return file_paths


def main():
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <input_file> <anomalies_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    anom_file = sys.argv[2]

    file_paths = generate_file_paths(input_file, anom_file)

    # Assuming script.py is the name of the script you want to call
    command = ["python3", "-m", "src.Analysis.plot_loader"] + file_paths
    os.system(" ".join(command))


if __name__ == "__main__":
    main()
