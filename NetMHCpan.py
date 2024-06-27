import sys
import pandas as pd
import mhctools
from mhctools import NetMHCpan
from mhcnames import normalize_allele_name

if len(sys.argv) < 4:
    print("Error: Missing arguments. Usage: mhctools_script.py <tool> <peptides_file> <alleles_file>")
    sys.exit(1)

tool = sys.argv[1]
peptides_file = sys.argv[2]
alleles_file = sys.argv[3]

print("Python version:", sys.version)
print("Pandas version:", pd.__version__)
print("mhctools version:", mhctools.__version__)

if tool == "NetMHCpan":
    try:
        # Read alleles from file
        with open(alleles_file, 'r') as f:
            alleles = [normalize_allele_name(line.strip()) for line in f if not line.startswith("#") and line.strip()]

        # Read peptides from file
        protein_sequences = {}
        with open(peptides_file, 'r') as f:
            for line in f:
                if not line.startswith("#") and line.strip():
                    parts = line.strip().split()
                    if len(parts) == 2:
                        protein_sequences[parts[0]] = parts[1]

        # Initialize NetMHCpan predictor
        pan_predictor = NetMHCpan(alleles=alleles)
        print("NetMHCpan loaded successfully.")

        # Predict binding
        binding_predictions = pan_predictor.predict_subsequences(protein_sequences, peptide_lengths=[9])

        # Convert predictions to DataFrame
        df = binding_predictions.to_dataframe()

        # Print strong binders
        for binding_prediction in binding_predictions:
            if binding_prediction.affinity < 100:
                print("Strong binder: %s" % (binding_prediction,))
    except Exception as e:
        print(f"Error running NetMHCpan: {e}")
else:
    print("Error: Invalid tool specified. Please specify 'NetMHCpan'.")
    sys.exit(1)

print("All modules imported successfully!")
