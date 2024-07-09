import subprocess
import pandas as pd

# Create a custom MHC sequence file
custom_mhc_file = "custom_mhc.fasta"
with open(custom_mhc_file, "w") as f:
    f.write(""">MHC1
MSAQRVGSLADGRTVEALHGAEGLRQSLPDC
>MHC2
MSLQRVGSLADGRTVEALHGAEGLRQSLPDC
""")

# Create a peptide sequence file
peptide_file = "peptides.fasta"
with open(peptide_file, "w") as f:
    f.write(""">1L2Y
NLYIQWLKDGGPSSGRPPPS
>1L3Y
ECDTINCERYNGQVCGGPGRGLCFCGKCRCHPGFEGSACQA
""")

# Run NetMHCpan with custom MHC sequences
result_file = "netmhcpan_custom_predictions.csv"
command = [
    "netMHCpan",
    "-f", peptide_file,
    "-inptype", "1",
    "-a", custom_mhc_file,
    "-xls",
    "-xlsfile", result_file
]

subprocess.run(command, check=True)

# Read and print the result
df = pd.read_csv(result_file, sep='\t')
print(df)

# Example: Print strong binders
strong_binders = df[df['Affinity(nM)'] < 100]
print("Strong binders:")
print(strong_binders)
