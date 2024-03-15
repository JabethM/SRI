import os
import sys
import plot_loader


def parse_filename(filename):
    """Parse numbers from filename."""
    parts = filename.split('_')[-1].split('.')[:3]
    return tuple(map(int, parts))


def get_anomalies(b, c, filename):
    """Extract existing numbers from the file after '- Anomalous files:'."""
    existing_numbers = []
    start_reading = False  # Set True to exclude single files
    with open(filename, 'r') as file:
        for line in file:
            if line.strip() == '- Single Files:':
                start_reading = True
                continue
            if line.strip() == '- Anomalous files:':
                start_reading = True  # Set False to include Anomalous files
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
    for i in repeats:
        version_parts[0] = str(i)
        new_filename = base_name.split('_')[0] + '_' + '.'.join(version_parts) + extension
        file_paths.append(os.path.join(directory, new_filename))

    return file_paths


def main():
    if len(sys.argv) != 4:
        print(f"Usage: python {sys.argv[0]} <plotData_file> <anomalies_file> <regime>")
        sys.exit(1)

    input_file = sys.argv[1]
    anom_file = sys.argv[2]
    regime = int(sys.argv[3])

    plot_paths = generate_file_paths(input_file, anom_file)
    if regime == 0:
        # Plot infectious curve.vs.time
        plot_loader.plot_PlotPickle_files(plot_paths)
    elif regime == 1:
        # Plot Connections.vs.time
        connection_paths = [p.split('plotData')[0] + 'connectionStatsData' + p.split('plotData')[1] for p in plot_paths]
        plot_loader.plot_ConnectionPickle_files(connection_paths)
    else:
        print("Regime OutOfBounds")
        sys.exit(1)


if __name__ == "__main__":
    main()
