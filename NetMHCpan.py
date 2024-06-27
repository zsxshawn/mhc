import sys
import pandas as pd
import mhctools
from mhctools import NetMHCpan, NetMHCIIpan
from mhcnames import normalize_allele_name

if len(sys.argv) < 2:
    print("Error: No tool specified. Please specify either 'NetMHCpan' or 'NetMHCIIpan'.")
    sys.exit(1)

tool = sys.argv[1]

print("Python version:", sys.version)
print("Pandas version:", pd.__version__)
print("mhctools version:", mhctools.__version__)

if tool == "NetMHCpan":
    try:
        pan_predictor = NetMHCpan(alleles=[normalize_allele_name("HLA-A*02:01")])
        print("NetMHCpan loaded successfully.")
    except Exception as e:
        print(f"Error loading NetMHCpan: {e}")
elif tool == "NetMHCIIpan":
    try:
        ii_pan_predictor = NetMHCIIpan(alleles=[normalize_allele_name("HLA-DRB1*01:01")])
        print("NetMHCIIpan loaded successfully.")
    except Exception as e:
        print(f"Error loading NetMHCIIpan: {e}")
else:
    print("Error: Invalid tool specified. Please specify either 'NetMHCpan' or 'NetMHCIIpan'.")
    sys.exit(1)

print("All modules imported successfully!")
