# Make a list of *.root files in a directory
# running example: python script.py /path/to/directory output.txt
# -- output.txt file will be created in ../data/sample_list/
import os
import argparse

def list_root_files(directory):
    root_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".root"):
                root_files.append(os.path.join(root, file))
    return root_files

def main():
    parser = argparse.ArgumentParser(description="List all .root files in a directory with absolute paths and save to a text file.")
    parser.add_argument("directory", type=str, help="Path to the directory to scan.")
    parser.add_argument("output", type=str, help="Path to the output text file.")
    args = parser.parse_args()
    
    directory = os.path.abspath(args.directory)
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        return
    
    root_files = list_root_files(directory)
    
    with open("../data/sample_list/" + args.output, "w") as f:
        if root_files:
            for file in root_files:
                f.write(file + "\n")
        else:
            f.write("No .root files found.\n")
    
    print(f"Results saved to {args.output}")

if __name__ == "__main__":
    main()
