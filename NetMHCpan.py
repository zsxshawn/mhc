import os
import sys
import subprocess
import pandas as pd
from datetime import datetime
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def read_alleles_and_sequences(file_path):
    logger.debug(f"Reading alleles from file: {file_path}")
    alleles_and_sequences = {}
    current_allele = None
    current_sequence = []
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_allele:
                        alleles_and_sequences[current_allele] = ''.join(current_sequence)
                    current_allele = line[1:]
                    current_sequence = []
                elif current_allele:
                    current_sequence.append(line)
        
        if current_allele:
            alleles_and_sequences[current_allele] = ''.join(current_sequence)
        
        logger.debug(f"Read {len(alleles_and_sequences)} alleles")
        logger.debug(f"Alleles: {list(alleles_and_sequences.keys())}")
    except Exception as e:
        logger.error(f"Error reading alleles file: {e}")
    return alleles_and_sequences

def read_peptides(file_path):
    logger.debug(f"Reading peptides from file: {file_path}")
    protein_sequences = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith("#") and line.strip():
                    parts = line.strip().split()
                    if len(parts) == 2:
                        protein_sequences[parts[0]] = parts[1]
        logger.debug(f"Read {len(protein_sequences)} protein sequences")
        logger.debug(f"Protein names: {list(protein_sequences.keys())}")
    except Exception as e:
        logger.error(f"Error reading peptides file: {e}")
    return protein_sequences

def run_netmhcpan(peptides_file, alleles_file, output_file):
    logger.debug("Running NetMHCpan")
    cmd = [
        "netMHCpan",
        "-p", peptides_file,
        "-a", "file://" + alleles_file,
        "-xls",
        "-xlsfile", output_file
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.debug("NetMHCpan completed successfully")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running NetMHCpan: {e}")
        logger.error(f"NetMHCpan stdout: {e.stdout}")
        logger.error(f"NetMHCpan stderr: {e.stderr}")
        raise

def predict_binding(tool, peptides_file, alleles_file, output_name=None):
    if tool == "NetMHCpan":
        try:
            output_folder = "output"
            os.makedirs(output_folder, exist_ok=True)
            
            if output_name:
                output_file = os.path.join(output_folder, f"{output_name}.xls")
            else:
                now = datetime.now().strftime("%Y%m%d_%H%M%S")
                output_file = os.path.join(output_folder, f"{now}_{tool}_predictions.xls")
            
            run_netmhcpan(peptides_file, alleles_file, output_file)
            
            logger.info(f"Predictions saved to {output_file}")
            
            # Read the output file and create a DataFrame
            df = pd.read_csv(output_file, sep='\t', skiprows=1)
            logger.debug(f"Created DataFrame with {len(df)} rows")
            
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
    
    result_df = predict_binding(tool, peptides_file, alleles_file, output_name)
    logger.info("Script completed")
