import pandas as pd

#reading excel file
file = pd.read_excel(r"D:\CHUM_git\16s_data\Microbiota_10_MiSeqReadSet_2020-12-08.xlsx")

#separating the different primers sequences
r1_primers = file["Séquence de l'amorce sens"][1].split(";")
r2_primers = file["Séquence de l'amorce antisens"][1].split(";")

#writing the primer sequences in separate file
# Specify the file path
fasta_file_path = r"D:\CHUM_git\16s_data\r1_primers.fasta"

# Open the file in write mode
with open(fasta_file_path, "w") as fasta_file:
    # Write each string to the file in FASTA format
    for idx, string in enumerate(r1_primers, start=1):
        fasta_file.write(f">primer{idx}"+"\n") #Write fasta header
        fasta_file.write(string + "\n")  # Write the sequence


# Specify the file path
fasta_file_path = r"D:\CHUM_git\16s_data\r2_primers.fasta"

# Open the file in write mode
with open(fasta_file_path, "w") as fasta_file:
    # Write each string to the file in FASTA format
    for idx, string in enumerate(r2_primers, start=1):
        fasta_file.write(f">primer{idx}"+"\n") #Write fasta header
        fasta_file.write(string + "\n")  # Write the sequence


