import sys

def copy_file(input_file, output_file):
    try:
        # Open the input file in read mode
        with open(input_file, 'r') as f_input:
            # Read the content of the input file
            file_content = f_input.read()

        # Open the output file in write mode
        with open(output_file, 'w') as f_output:
            # Write the content to the output file
            f_output.write(file_content)

        print(f"Content copied from {input_file} to {output_file} successfully.")

    except FileNotFoundError:
        print("File not found. Please check the file path and try again.")
    except IOError:
        print("An error occurred while reading/writing the file.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Example usage:
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

copy_file(input_file_path, output_file_path)
