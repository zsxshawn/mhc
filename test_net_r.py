import subprocess
import pandas as pd
import tempfile
import os

def run_netmhcpan(allele, peptides, output_file):
    # Create a temporary file for peptides
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        for peptide in peptides:
            temp_file.write(f"{peptide}\n")
        temp_file_name = temp_file.name

    # Run NetMHCpan
    cmd = f"netMHCpan -a {allele} -p {temp_file_name} -xls -xlsfile {output_file}"
    subprocess.run(cmd, shell=True, check=True)

    # Remove temporary file
    os.unlink(temp_file_name)

    # Read and process the output
    df = pd.read_csv(output_file, sep='\t', skiprows=1)
    return df

# Example usage
allele = "HLA-A02:01"  # You can change this to any allele
peptides = [
    "NLYIQWLKD",
    "LYIQWLKDG",
    "YIQWLKDGG",
    "IQWLKDGGP",
    "QWLKDGGPS",
    "WLKDGGPSS",
    "LKDGGPSSG",
    "KDGGPSSGR",
    "DGGPSSGR"
]

output_file = "netmhcpan_results.xls"

results = run_netmhcpan(allele, peptides, output_file)

# Print strong binders
strong_binders = results[results['nM'] < 100]
print("Strong binders:")
print(strong_binders[['Peptide', 'nM', '%Rank']])

print(f"\nAll results saved to {output_file}")
