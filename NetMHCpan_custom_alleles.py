from mhctools import NetMHCpan
import pandas as pd

# New custom MHC sequence
custom_mhc_sequences = """>B*070201
MLVMAPRTVLLLLSAALALTETWAGSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRF
DSDAASPREEPRAPWIEQEGPEYWDRNTQIYKAQAQTDRESLRNLRGYYNQSEAGSHTLQ
SMYGCDVGPDGRLLRGHDQYAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQR
RAYLEGECVEWLRRYLENGKDKLERADPPKTHVTHHPISDHEATLRCWALGFYPAEITLT
WQRDGEDQTQDTELVETRPAGDRTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEP
SSQSTVPIVGIVAGLAVLAVVVIGAVVAAVMCRRKSSGGKGGSYSQAACSDSAQGSDVSL
TA"""

# Debug: Print custom MHC sequences to ensure they are correctly set
print("Custom MHC Sequences:\n", custom_mhc_sequences)

# Use only the custom MHC sequence
predictor = NetMHCpan(
    alleles=["Custom_B*070201"],
    custom_mhc_sequences=custom_mhc_sequences
)

# Scan the short proteins 1L2Y and 1L3Y for epitopes
protein_sequences = {
    "1L2Y": "NLYIQWLKDGGPSSGRPPPS",
    "1L3Y": "ECDTINCERYNGQVCGGPGRGLCFCGKCRCHPGFEGSACQA"
}

# Predict binding for subsequences of specified lengths
binding_predictions = predictor.predict_subsequences(protein_sequences, peptide_lengths=[9])

# Flatten binding predictions into a Pandas DataFrame
df = binding_predictions.to_dataframe()

# Write the DataFrame to a CSV file
df.to_csv("netmhcpan_custom_predictions.csv", index=False)

# Print the DataFrame to the console
print(df)

# Example: Print strong binders
for binding_prediction in binding_predictions:
    if binding_prediction.affinity < 100:
        print("Strong binder: %s" % (binding_prediction,))
