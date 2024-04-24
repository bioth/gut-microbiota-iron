import pandas as pd
file = pd.read_excel(r"D:\CHUM_git\16s_data\Microbiota_10_MiSeqReadSet_2020-12-08.xlsx")
print(file["Nom"])
print(len(pd.unique(file["Nom"])))
print(pd.duplicated(file))