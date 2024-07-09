from mhctools import NetMHCpan
import pandas as pd

# Custom MHC sequence
custom_mhc_sequences = """>USER_DEF
MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRF
DSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQ
RMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQL
RAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLT
WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEP
SSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSL
TACKV"""

# Protein sequences to be scanned
protein_sequences = {
    "PROTEIN1": "MKVLWAALLVTFLAGCQA",
    "PROTEIN2": "RWMVLWAALTVTLAGCQAI",
    "PROTEIN3": "KMVLWAALLVTFLAGCQAT"
}

# Create the predictor with the custom MHC sequence
predictor = NetMHCpan(
    alleles=["USER_DEF"],
    custom_mhc_sequences=custom_mhc_sequences
)

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
