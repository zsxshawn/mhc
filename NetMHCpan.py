import sys
import pandas as pd
import mhctools
from mhctools import NetMHCpan, NetMHCIIpan
from mhcnames import normalize_allele_name

print("Python version:", sys.version)
print("Pandas version:", pd.__version__)
print("mhctools version:", mhctools.__version__)

# Test creation of NetMHCpan and NetMHCIIpan predictors to check versions
try:
    pan_predictor = NetMHCpan(alleles=[normalize_allele_name("HLA-A*02:01")])
    print("NetMHCpan loaded successfully.")
except Exception as e:
    print(f"Error loading NetMHCpan: {e}")

try:
    ii_pan_predictor = NetMHCIIpan(alleles=[normalize_allele_name("HLA-DRB1*01:01")])
    print("NetMHCIIpan loaded successfully.")
except Exception as e:
    print(f"Error loading NetMHCIIpan: {e}")

print("All modules imported successfully!")
print("From Github")
