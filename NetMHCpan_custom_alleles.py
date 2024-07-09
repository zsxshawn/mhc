from mhctools import NetMHCpan
import pandas as pd

custom_mhc_sequences = """>MHC1
MSAQRVGSLADGRTVEALHGAEGLRQSLPDC
>MHC2
MSLQRVGSLADGRTVEALHGAEGLRQSLPDC
"""

# Run NetMHCpan for alleles HLA-A*01:01 and HLA-A*02:01 with custom MHC sequences
predictor = NetMHCpan(
    alleles=["A*02:01", "hla-a0101"],
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
