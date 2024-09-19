import pandas as pd
import sys
import os 

#Moving to metadata folder
folder_path = sys.argv[1]
directory = os.path.dirname(folder_path)
folder_path = os.path.join(directory,"metadata").replace("\\","/")



# Open the file in write mode
def writeAdaptersFasta(fasta_file_path, primers):
    with open(fasta_file_path, "w") as fasta_file:
        # Write each string to the file in FASTA format
        for idx, string in enumerate(primers, start=1):
            fasta_file.write(f">primer{idx}"+"\n") #Write fasta header
            fasta_file.write(string + "\n")  # Write the sequence


def file_contains_string(filename, search_string):
    with open(filename, 'r') as file:
        for line in file:
            if search_string in line:
                return True
    return False


def find_file_containing_string(directory, search_string):
    for root, _, files in os.walk(directory):
        for filename in files:
            filepath = os.path.join(root, filename)
            if file_contains_string(filepath, search_string):
                return filepath
    return None

def main():
    print(directory)
    print(folder_path)

    #reading excel file
    file = pd.read_excel(find_file_containing_string(directory, "SeqReadSet"))

    #separating the different primers sequences
    r1_primers = file["Séquence de l'amorce sens"][1].split(";")
    r2_primers = file["Séquence de l'amorce antisens"][1].split(";")

    #writing the primer sequences in separate file
    # Specify the file path
    fasta_file_path1 = os.path.join(directory,"primers/r1_primers.fasta")

    fasta_file_path2 = os.path.join(directory,"primers/r2_primers.fasta")

    writeAdaptersFasta(fasta_file_path1, r1_primers)
    writeAdaptersFasta(fasta_file_path2, r2_primers)

# Call the main function if this script is run directly
if __name__ == "__main__":
    main()

