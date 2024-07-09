from mhctools import NetMHCpan4
import os

# Create a custom MHC sequence file
custom_mhc_file = "custom_mhc.fasta"
with open(custom_mhc_file, "w") as f:
    f.write(""">MHC1
MSAQRVGSLADGRTVEALHGAEGLRQSLPDC
>MHC2
MSLQRVGSLADGRTVEALHGAEGLRQSLPDC
""")

# Initialize NetMHCpan predictor with the custom MHC file
# Note: This step might not work if the library does not support custom MHC sequences directly
predictor = NetMHCpan4(alleles=[custom_mhc_file])

# Define your protein sequences
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
