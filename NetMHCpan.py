import os
import sys
import pandas as pd
import mhctools
from mhctools import NetMHCpan
from mhcnames import normalize_allele_name

def read_alleles_and_pseudo_sequences(file_path):
    alleles_and_sequences = {}
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                parts = line.strip().split()
                allele = normalize_allele_name(parts[0])
                if len(parts) == 2:
                    pseudo_sequence = parts[1]
                    alleles_and_sequences[allele] = pseudo_sequence
                else:
                    alleles_and_sequences[allele] = None
    return alleles_and_sequences

def read_peptides(file_path):
    protein_sequences = {}
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                parts = line.strip().split()
                if len(parts) == 2:
                    protein_sequences[parts[0]] = parts[1]
    return protein_sequences

def predict_binding(tool, peptides_file, alleles_file):
    if tool == "NetMHCpan":
        try:
            # Read alleles and pseudo sequences from file
            alleles_and_sequences = read_alleles_and_pseudo_sequences(alleles_file)
            
            # Read peptides from file
            protein_sequences = read_peptides(peptides_file)
            
            # Initialize NetMHCpan predictor with alleles
            pan_predictor = NetMHCpan(alleles=list(alleles_and_sequences.keys()))
            print("NetMHCpan loaded successfully.")
            
            # Predict binding
            binding_predictions = pan_predictor.predict_subsequences(protein_sequences, peptide_lengths=[9])
            
            # Convert predictions to DataFrame
            df = binding_predictions.to_dataframe()
            
            # Print strong binders
            for binding_prediction in binding_predictions:
                if binding_prediction.affinity < 100:
                    print(f"Strong binder: {binding_prediction}")
                    if alleles_and_sequences[binding_prediction.allele]:
                        print(f"Pseudo sequence: {alleles_and_sequences[binding_prediction.allele]}")
            
            return df
        
        except Exception as e:
            print(f"Error running NetMHCpan: {e}")
    else:
        print("Error: Invalid tool specified. Please specify 'NetMHCpan'.")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Error: Missing arguments. Usage: NetMHCpan.py <tool> <peptides_file> <alleles_file>")
        sys.exit(1)
    
    tool = sys.argv[1]
    peptides_file = sys.argv[2]
    alleles_file = sys.argv[3]
    
    print("Python version:", sys.version)
    print("mhctools version:", mhctools.__version__)
    
    result_df = predict_binding(tool, peptides_file, alleles_file)
    print(result_df)
