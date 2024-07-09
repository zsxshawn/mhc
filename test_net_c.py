import subprocess
import pandas as pd
import tempfile
import os

def run_netmhcpan(allele, peptides, output_file, mhc_sequence=None):
    # Create a temporary file for peptides
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        for peptide in peptides:
            temp_file.write(f"{peptide}\n")
        temp_file_name = temp_file.name

    # Base command
    cmd = f"netMHCpan -p {temp_file_name} -xls -xlsfile {output_file}"

    # If MHC sequence is provided, create a temporary file for it
    if mhc_sequence:
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as mhc_file:
            mhc_file.write(f">Custom_MHC\n{mhc_sequence}\n")
            mhc_file_name = mhc_file.name
        cmd += f" -inptype 1 -a {mhc_file_name}"
    else:
        cmd += f" -a {allele}"

    # Run NetMHCpan
    subprocess.run(cmd, shell=True, check=True)

    # Remove temporary files
    os.unlink(temp_file_name)
    if mhc_sequence:
        os.unlink(mhc_file_name)

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

# For custom MHC sequence
custom_mhc = """MLVMAPRTVLLLLSAALALTETWAGSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRF
DSDAASPREEPRAPWIEQEGPEYWDRNTQIYKAQAQTDRESLRNLRGYYNQSEAGSHTLQ
SMYGCDVGPDGRLLRGHDQYAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQR
RAYLEGECVEWLRRYLENGKDKLERADPPKTHVTHHPISDHEATLRCWALGFYPAEITLT
WQRDGEDQTQDTELVETRPAGDRTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEP
SSQSTVPIVGIVAGLAVLAVVVIGAVVAAVMCRRKSSGGKGGSYSQAACSDSAQGSDVSL
TA"""

# Choose whether to use allele name or custom MHC sequence
use_custom_mhc = True  # Set to False to use allele name instead

if use_custom_mhc:
    results = run_netmhcpan(None, peptides, output_file, mhc_sequence=custom_mhc)
else:
    results = run_netmhcpan(allele, peptides, output_file)

# Print all results
print("All results:")
print(results)

# Print strong binders (using EL_Rank < 0.5 as the threshold for strong binders)
strong_binders = results[results['EL_Rank'] < 0.5]
print("\nStrong binders:")
if not strong_binders.empty:
    print(strong_binders[['Peptide', 'EL-score', 'EL_Rank']])
else:
    print("No strong binders found.")

print(f"\nAll results saved to {output_file}")
