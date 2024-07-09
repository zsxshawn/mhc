from mhctools import NetMHCpan
import pandas as pd

custom_mhc_sequences = """>MHC1
MSAQRVGSLADGRTVEALHGAEGLRQSLPDC
>MHC2
MSLQRVGSLADGRTVEALHGAEGLRQSLPDC
"""

# Run NetMHCpan with custom MHC sequences
predictor = NetMHCpan(custom_mhc_sequences=custom_mhc_sequences)

# Define protein sequences
protein_sequences = {
    "1L2Y": "NLYIQWLKDGGPSSGRPPPS",
    "1L3Y": "ECDTINCERYNGQVCGGPGRGLCFCGKCRCHPGFEGSACQA"
}

# Predict binding for subsequences of specified lengths
binding_predictions = predictor.predict_subsequences(protein_sequences, peptide_lengths=[9])

# Flatten binding predictions into a Pandas DataFrame
df = binding_predictions.to_dataframe()

# Save the DataFrame to a CSV file
result_file = 'netmhcpan_custom_predictions.csv'
df.to_csv(result_file, sep='\t', index=False)

# Print the DataFrame
print(df)
