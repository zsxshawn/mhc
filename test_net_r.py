import subprocess
import pandas as pd
import tempfile
import os

def run_netmhcpan(peptides, output_file, mhc_sequence):
    # Create a temporary file for peptides
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        for peptide in peptides:
            temp_file.write(f"{peptide}\n")
        peptide_file_name = temp_file.name

    # Prepare MHC sequence (remove newlines and spaces)
    mhc_sequence_cleaned = ''.join(mhc_sequence.split())

    # Run NetMHCpan
    cmd = f"netMHCpan -p {peptide_file_name} -xls -xlsfile {output_file} -sequence {mhc_sequence_cleaned}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    # Print the command and its output for debugging
    print(f"Command: {cmd}")
    print("Standard Output:")
    print(result.stdout)
    print("Standard Error:")
    print(result.stderr)

    # Remove temporary file
    os.unlink(peptide_file_name)

    # Read and process the output
    if os.path.exists(output_file):
        df = pd.read_csv(output_file, sep='\t', skiprows=1)
        return df
    else:
        print(f"Error: Output file {output_file} was not created.")
        return None

# Example usage
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

# Custom MHC sequence (B*070201)
custom_mhc = """MLVMAPRTVLLLLSAALALTETWAGSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRF
DSDAASPREEPRAPWIEQEGPEYWDRNTQIYKAQAQTDRESLRNLRGYYNQSEAGSHTLQ
SMYGCDVGPDGRLLRGHDQYAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQR
RAYLEGECVEWLRRYLENGKDKLERADPPKTHVTHHPISDHEATLRCWALGFYPAEITLT
WQRDGEDQTQDTELVETRPAGDRTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEP
SSQSTVPIVGIVAGLAVLAVVVIGAVVAAVMCRRKSSGGKGGSYSQAACSDSAQGSDVSL
TA"""
results = run_netmhcpan(peptides, output_file, mhc_sequence=custom_mhc)

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
