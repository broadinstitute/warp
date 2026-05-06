

META_DICT_BASE = {
    "filter": {
        "artifact_prone_site": {
            "Description": "Variant overlaps an artifact-prone site"
        }
    },
    "format": {
        "DP": {"Description": "Depth of coverage", "Number": "1", "Type": "Integer"},
        "FT": {
            "Description": "Sample-level genotype filters",
            "Number": ".",
            "Type": "String",
        },
        "HL": {"Description": "Heteroplasmy level", "Number": "1", "Type": "Float"},
        "MQ": {"Description": "Mapping quality", "Number": "1", "Type": "Float"},
        "TLOD": {
            "Description": "Log 10 likelihood ratio score of variant existing versus not existing",
            "Number": "1",
            "Type": "Float",
        },
    },
}

META_DICT_V2_FMT = {
    'AD': {"Description": "Allelic depth of REF and ALT", "Number": "R", "Type": "Integer"},
    'OriginalSelfRefAlleles': {
        'Description':'Original self-reference alleles (only if alleles were changed in Liftover repair pipeline)', 
        'Number':'R', 
        'Type':'String'},
    'SwappedFieldIDs': {
        'Description':'Fields remapped during liftover (only if alleles were changed in Liftover repair pipeline)', 
        'Number':'1', 
        'Type':'String'
    },
    'F2R1': {"Description": "Count of reads in F2R1 pair orientation supporting each allele", "Number": "R", "Type": "Integer"},
    'F1R2': {"Description": "Count of reads in F1R2 pair orientation supporting each allele", "Number": "R", "Type": "Integer"},
    "MPOS": {
        "Description": "median distance from end of read",
        "Number": "1",
        "Type": "Integer",
    },
    "AS_SB_TABLE": {
        "Description": "Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.",
        "Number": "1",
        "Type": "String",
    },
    "STR": {
        "Description": "Variant is a short tandem repeat. 1 if True, 0 if False.",
        "Number": "1",
        "Type": "Integer",
    },
    "STRQ": {
        "Description": "Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors",
        "Number": "1",
        "Type": "Integer",
    },
    "RPA": {
        "Description": "Number of times tandem repeat unit is repeated, for each allele (including reference)",
        "Number": "R",
        "Type": "Integer",
    },
}

V2_FIELD_KEY = {
    'OriginalSelfRefAlleles':'array<str>', 
    'SwappedFieldIDs':'str',
    'F2R1':'array<int32>', 
    'F1R2':'array<int32>',
    'MPOS':'int32',
    'AS_SB_TABLE':'str',
    'STR':'int32',
    'STRQ':'int32',
    'RPA':'array<int32>'
}

V2_INFO_TO_FORMAT = ['MPOS','AS_SB_TABLE','STR','STRQ','RPA']
V2_INFO_TO_FORMAT_REQUIREINDEX = ['MPOS']


LIFTOVERFILTERS = set(['NoTarget','MismatchedRefAllele','IndelStraddlesMultipleIntevals'])
CUSTOMLIFTOVERFILTERS = set(['FailedPicardLiftoverVcf', 'InsertionSharesForceCalledDeletion', 'InsertionSharesForceCalledInsertion',
                             'AddGRCh38RefDeleToRefSiteIns', 'ComplexSwapField', 'NewInsertionHaplotype',
                             'SwapFirstAlleleIndel', 'ReplaceInternalBaseDeletion', 'FancyFieldInversion',
                             'DeletionSpannedHomoplasmicInsertion', 'LiftoverSuccessEntrySwap', 'ForceCalledHomoplasmy',
                             'LeftShiftedIndel', 'FailedDuplicateVariant'])