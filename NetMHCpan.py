import os
import sys
import pandas as pd
from datetime import datetime
from mhctools import NetMHCpan41
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def read_alleles_and_pseudo_sequences(file_path):
    logger.info(f"Reading alleles from file: {file_path}")
    alleles_and_sequences = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith("#") and line.strip():
                    parts = line.strip().split()
                    allele = parts[0]
                    if len(parts) == 2:
                        pseudo_sequence = parts[1]
                        alleles_and_sequences[allele] = pseudo_sequence
                    else:
                        alleles_and_sequences[allele] = None
        logger.info(f"Read {len(alleles_and_sequences)} alleles")
        logger.info(f"Alleles: {list(alleles_and_sequences.keys())}")
    except Exception as e:
        logger.error(f"Error reading alleles file: {e}")
    return alleles_and_sequences

def read_peptides(file_path):
    logger.info(f"Reading peptides from file: {file_path}")
    protein_sequences = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith("#") and line.strip():
                    parts = line.strip().split()
                    if len(parts) == 2:
                        protein_sequences[parts[0]] = parts[1]
        logger.info(f"Read {len(protein_sequences)} protein sequences")
        logger.info(f"Protein names: {list(protein_sequences.keys())}")
    except Exception as e:
        logger.error(f"Error reading peptides file: {e}")
    return protein_sequences

def predict_binding(tool, peptides_file, alleles_file, output_name=None):
    logger.info(f"Starting prediction with {tool}")
    if tool == "NetMHCpan":
        try:
            alleles_and_sequences = read_alleles_and_pseudo_sequences(alleles_file)
            protein_sequences = read_peptides(peptides_file)
            
            logger.info("Initializing NetMHCpan predictor")
            predictor = NetMHCpan41(
                alleles=list(alleles_and_sequences.keys()),
                program_name="netMHCpan",
                process_limit=-1,
                default_peptide_lengths=[9],
                extra_flags=["-p"] + [f"{allele}:{seq}" for allele, seq in alleles_and_sequences.items() if seq]
            )
            logger.info("NetMHCpan predictor initialized successfully")
            
            logger.info("Starting prediction")
            binding_predictions = predictor.predict_subsequences(protein_sequences)
            logger.info(f"Prediction complete. Found {len(binding_predictions)} predictions")
            
            df = binding_predictions.to_dataframe()
            logger.info(f"Created DataFrame with {len(df)} rows")
            
            for binding_prediction in binding_predictions:
                if binding_prediction.affinity < 100:
                    logger.info(f"Strong binder: {binding_prediction}")
                    if alleles_and_sequences[binding_prediction.allele]:
                        logger.info(f"Pseudo sequence: {alleles_and_sequences[binding_prediction.allele]}")
            
            output_folder = "output"
            os.makedirs(output_folder, exist_ok=True)
            
            if output_name:
                output_file = os.path.join(output_folder, f"{output_name}.csv")
            else:
                now = datetime.now().strftime("%Y%m%d_%H%M%S")
                output_file = os.path.join(output_folder, f"{now}_{tool}_predictions.csv")
            
            df.to_csv(output_file, index=False)
            logger.info(f"Predictions saved to {output_file}")
            
            return df
        
        except Exception as e:
            logger.error(f"Error running NetMHCpan: {e}")
            logger.exception("Exception details:")
    else:
        logger.error("Error: Invalid tool specified. Please specify 'NetMHCpan'.")
        sys.exit(1)

if __name__ == "__main__":
    logger.info("Script started")
    if len(sys.argv) < 4:
        logger.error("Error: Missing arguments. Usage: NetMHCpan.py <tool> <peptides_file> <alleles_file> [output_name]")
        sys.exit(1)
    
    tool = sys.argv[1]
    peptides_file = sys.argv[2]
    alleles_file = sys.argv[3]
    output_name = sys.argv[4] if len(sys.argv) > 4 else None
    
    logger.info(f"Arguments: tool={tool}, peptides_file={peptides_file}, alleles_file={alleles_file}, output_name={output_name}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"mhctools version: {mhctools.__version__}")
    
    result_df = predict_binding(tool, peptides_file, alleles_file, output_name)
    logger.info("Script completed")
