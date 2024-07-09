import os
import sys
import subprocess
import pandas as pd
from datetime import datetime
import logging
import shutil

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def check_netmhcpan_installation():
    netmhcpan_path = shutil.which("netMHCpan")
    if netmhcpan_path:
        logger.info(f"NetMHCpan found at: {netmhcpan_path}")
    else:
        logger.error("NetMHCpan not found in PATH")

def run_netmhcpan(peptides_file, alleles_file, output_file):
    logger.debug("Running NetMHCpan")
    cmd = [
        "netMHCpan",
        "-p", os.path.abspath(peptides_file),
        "-a", f"file://{os.path.abspath(alleles_file)}",
        "-xls",
        "-xlsfile", os.path.abspath(output_file)
    ]
    
    logger.debug(f"Executing command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.debug("NetMHCpan completed successfully")
        logger.debug(f"NetMHCpan stdout: {result.stdout}")
        if result.stderr:
            logger.warning(f"NetMHCpan stderr: {result.stderr}")
        
        # Check if the output file exists
        if os.path.exists(output_file):
            logger.debug(f"Output file created successfully: {output_file}")
        else:
            logger.error(f"Output file not found: {output_file}")
        
        # Print contents of the current directory and output directory
        logger.debug(f"Contents of current directory:")
        for file in os.listdir():
            logger.debug(f"  {file}")
        
        output_dir = os.path.dirname(output_file)
        logger.debug(f"Contents of output directory {output_dir}:")
        for file in os.listdir(output_dir):
            logger.debug(f"  {file}")
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running NetMHCpan: {e}")
        logger.error(f"NetMHCpan stdout: {e.stdout}")
        logger.error(f"NetMHCpan stderr: {e.stderr}")
        raise

def predict_binding(tool, peptides_file, alleles_file, output_name=None):
    if tool == "NetMHCpan":
        try:
            output_folder = os.path.abspath("output")
            os.makedirs(output_folder, exist_ok=True)
            
            if output_name:
                output_file = os.path.join(output_folder, f"{output_name}.xls")
            else:
                now = datetime.now().strftime("%Y%m%d_%H%M%S")
                output_file = os.path.join(output_folder, f"{now}_{tool}_predictions.xls")
            
            logger.debug(f"Output file will be: {output_file}")
            
            # Check write permissions
            if os.access(output_folder, os.W_OK):
                logger.debug(f"Have write permissions for {output_folder}")
            else:
                logger.error(f"No write permissions for {output_folder}")
            
            run_netmhcpan(peptides_file, alleles_file, output_file)
            
            if os.path.exists(output_file):
                logger.info(f"Predictions saved to {output_file}")
                
                # Read the output file and create a DataFrame
                df = pd.read_csv(output_file, sep='\t', skiprows=1)
                logger.debug(f"Created DataFrame with {len(df)} rows")
                
                return df
            else:
                logger.error(f"Output file not created: {output_file}")
                return None
        
        except Exception as e:
            logger.error(f"Error running NetMHCpan: {e}")
            logger.exception("Exception details:")
            return None
    else:
        logger.error("Error: Invalid tool specified. Please specify 'NetMHCpan'.")
        sys.exit(1)

if __name__ == "__main__":
    logger.info("Script started")
    logger.info(f"Current working directory: {os.getcwd()}")
    check_netmhcpan_installation()
    
    if len(sys.argv) < 4:
        logger.error("Error: Missing arguments. Usage: NetMHCpan.py <tool> <peptides_file> <alleles_file> [output_name]")
        sys.exit(1)
    
    tool = sys.argv[1]
    peptides_file = os.path.abspath(sys.argv[2])
    alleles_file = os.path.abspath(sys.argv[3])
    output_name = sys.argv[4] if len(sys.argv) > 4 else None
    
    logger.info(f"Arguments: tool={tool}, peptides_file={peptides_file}, alleles_file={alleles_file}, output_name={output_name}")
    logger.info(f"Python version: {sys.version}")
    
    result_df = predict_binding(tool, peptides_file, alleles_file, output_name)
    logger.info("Script completed")
