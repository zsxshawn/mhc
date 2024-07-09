from mhctools import NetMHCpan

DEFAULT_ALLELE = 'HLA-A*02:01'

protein_sequence_dict = {
    "SMAD4-001": "ASIINFKELA",
    "TP53-001": "ASILLLVFYW"
}

# Tests will also be run using the following program names if they are installed.
# In any case, a program called "netMHCpan" MUST be installed and working for
# this test suite to succeed.
OPTIONAL_NETMHCPAN_PROGRAM_NAMES = [
    "netMHCpan-2.8",
    "netMHCpan-3.0",
    "netMHCpan-4.0",
    "netMHCpan-4.1",
]

# Custom MHC sequences in FASTA format
custom_mhc_sequences = """>MHC1
MSAQRVGSLADGRTVEALHGAEGLRQSLPDC
>MHC2
MSLQRVGSLADGRTVEALHGAEGLRQSLPDC
"""

def test_netmhc_pan():
    check_netmhc_pan("netMHCpan", True)  # required

    for program_name in OPTIONAL_NETMHCPAN_PROGRAM_NAMES:
        check_netmhc_pan(program_name, False)  # optional

def check_netmhc_pan(program_name, fail_if_no_such_program=True):
    try:
        predictor = NetMHCpan(
            alleles=[DEFAULT_ALLELE], program_name=program_name)
    except FileNotFoundError:
        if fail_if_no_such_program:
            raise
        print("Skipping because no such program: %s" % program_name)
        return

    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[9])
    assert len(binding_predictions) == 4, \
        "Expected 4 binding predictions from %s" % (binding_predictions,)
    for x in binding_predictions:
        offset = x.offset
        length = x.length
        seq_name = x.source_sequence_name
        expected_peptide = protein_sequence_dict[seq_name][offset:offset + length]
        eq_(expected_peptide, x.peptide,
            "Peptide mismatch: expected %s but got %s in binding prediction '%s'" % (
                expected_peptide, x.peptide, x,))

def test_netmhc_pan_custom_mhc():
    predictor = NetMHCpan(custom_mhc_sequences=custom_mhc_sequences, program_name="netMHCpan-4.1")

    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[9])
    assert len(binding_predictions) == 4, \
        "Expected 4 binding predictions from %s" % (binding_predictions,)
    for x in binding_predictions:
        offset = x.offset
        length = x.length
        seq_name = x.source_sequence_name
        expected_peptide = protein_sequence_dict[seq_name][offset:offset + length]
        eq_(expected_peptide, x.peptide,
            "Peptide mismatch: expected %s but got %s in binding prediction '%s'" % (
                expected_peptide, x.peptide, x,))

def test_netmhc_pan_multiple_lengths():
    predictor = NetMHCpan(alleles=["A6801"])
    binding_predictions = predictor.predict_peptides(
        ["A" * 8, "A" * 9, "A" * 10, "A" * 11])
    assert len(binding_predictions) == 4, \
        "Expected 4 epitopes from %s" % (binding_predictions,)

def test_netmhc_pan_multiple_alleles():
    alleles = 'A*02:01,B*35:02'
    predictor = NetMHCpan(
        alleles=alleles,
        default_peptide_lengths=[9])
    sequence_dict = {
        "SMAD4-001": "ASIINFKELA",
        "TP53-001": "ASILLLVFYW"
    }
    binding_predictions = predictor.predict_subsequences(
        sequence_dict=sequence_dict)
    assert len(binding_predictions) == 8, \
        "Expected 8 binding predictions from %s" % (binding_predictions,)
