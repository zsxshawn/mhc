from mhctools import NetMHCpan

# Custom MHC sequences in FASTA format
custom_mhc_sequences = """>MHC1
MSAQRVGSLADGRTVEALHGAEGLRQSLPDC
>MHC2
MSLQRVGSLADGRTVEALHGAEGLRQSLPDC
"""

# Initialize NetMHCpan predictor with custom MHC sequences
predictor = NetMHCpan(custom_mhc_sequences=custom_mhc_sequences)

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
