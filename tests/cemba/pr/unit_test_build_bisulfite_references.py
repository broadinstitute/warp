#!/usr/bin/env python2

import filecmp
import os
from scripts.build_bisulfite_references import build_ref
import unittest
import tempfile


class BuildRefTest(unittest.TestCase):

    def test_files(self):
        # file paths for truth files
        current_file_dir = os.path.dirname(__file__)
        truth_input_ref_fasta = os.path.join(current_file_dir, "truth_input_ref.fa")
        truth_fwd_convert_ref_fasta = os.path.join(current_file_dir, "truth_fwd_converted_ref.fa")
        truth_rev_convert_ref_fasta =  os.path.join(current_file_dir, "truth_rev_converted_ref.fa")

        # generate temp file names for output files from build references python script
        fwd_converted_reference_fasta = os.path.join(tempfile.mkdtemp(), 'fwd_converted_reference_fasta')
        rev_converted_reference_fasta = os.path.join(tempfile.mkdtemp(), 'rev_converted_reference_fasta')

        # run the build references python script with given ref fasta
        build_ref(truth_input_ref_fasta, fwd_converted_reference_fasta, rev_converted_reference_fasta)

        # compare expected outputs from python script to truth files
        self.assertTrue(filecmp.cmp(fwd_converted_reference_fasta,
                                    truth_fwd_convert_ref_fasta,
                                    shallow=False), "Error: The forward converted file is not equal to truth file")
        self.assertTrue(filecmp.cmp(rev_converted_reference_fasta,
                                    truth_rev_convert_ref_fasta,
                                    shallow=False), "Error: The reverse converted file is not equal to truth file")


if __name__ == '__main__':
    unittest.main()
