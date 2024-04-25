import pandas as pd

# Open the file in write mode
def writeAdaptersFasta(fasta_file_path, primers):
    with open(fasta_file_path, "w") as fasta_file:
        # Write each string to the file in FASTA format
        for idx, string in enumerate(primers, start=1):
            fasta_file.write(f">primer{idx}"+"\n") #Write fasta header
            fasta_file.write(string + "\n")  # Write the sequence

def main():
    writeAdaptersFasta(fasta_file_path1, r1_primers)
    writeAdaptersFasta(fasta_file_path2, r2_primers)
    #reading excel file
    file = pd.read_excel(r"../../../16s_data/Microbiota_10_MiSeqReadSet_2020-12-08.xlsx")

    #separating the different primers sequences
    r1_primers = file["Séquence de l'amorce sens"][1].split(";")
    r2_primers = file["Séquence de l'amorce antisens"][1].split(";")

    #writing the primer sequences in separate file
    # Specify the file path
    fasta_file_path1 = r"../../../16s_data/r1_primers.fasta"
    fasta_file_path2 = r"../../../16s_data/r2_primers.fasta"

# Call the main function if this script is run directly
if __name__ == "__main__":
    main()

