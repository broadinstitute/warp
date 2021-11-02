#!/usr/bin/env python

import io
import argparse


def build_ref(input_fasta, fwd_filename, rev_filename):
    """
    builds the bisulfite references to use with bowtie2
    input_fasta: the reference fasta to be converted
    fwd_filename: the filename/path for the converted forward fasta (c's to t's)
    rev_filename: the filename for the converted backward/reverse fasta (g's to a's)
    """

    # convert the bases in from the input fasta
    with open(fwd_filename, "w") as fwd_out, open(rev_filename, "w") as rev_out:
        with open(input_fasta, "r") as ref_in:
            # read each character/base via buffer instead of by line (new line character is unreliable)
            in_header = False
            write_header_done = False
            buffer = ref_in.read(io.DEFAULT_BUFFER_SIZE)
            while buffer:
                for base in buffer:
                    if in_header:
                        if base == " " and not write_header_done:
                            fwd_out.write("_CT_converted\n")
                            rev_out.write("_GA_converted\n")
                            write_header_done = True

                        if base == "\n":
                            if not write_header_done:
                                fwd_out.write("_CT_converted\n")
                                rev_out.write("_GA_converted\n")
                                write_header_done = True

                            in_header = False

                        if not write_header_done:
                            fwd_out.write(base)
                            rev_out.write(base)

                    else:
                        fwd_out.write(get_converted_base(base, "c", "T"))
                        rev_out.write(get_converted_base(base, "g", "A"))

                        if base == ">":
                            in_header = True
                            write_header_done = False

                buffer = ref_in.read(io.DEFAULT_BUFFER_SIZE)


def get_converted_base(base, unconverted_base, converted_base):
    """
    :param base: base that is subject to conversion
    :param unconverted_base: base that should be converted
    :param converted_base: base that should be converted to
    """

    #  if the base is a c/C convert it to a T
    if base.upper() == unconverted_base.upper():
        return converted_base.upper()

    # capitalize base
    else:
        return base.upper()


if __name__ == "__main__":
    # get the argument inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fasta",
                        "-i",
                        dest="input_fasta",
                        required=True,
                        help="the input reference fasta")
    parser.add_argument("--forward-convert-out",
                        "-f",
                        dest="forward_convert_out",
                        required=True,
                        help="the filename/path for the converted forward fasta (c's to t's)")
    parser.add_argument("--reverse-convert-out",
                        "-r",
                        dest="reverse_convert_out",
                        required=True,
                        help="the filename for the converted backward/reverse fasta (g's to a's)")
    args = parser.parse_args()

    build_ref(args.input_fasta, args.forward_convert_out, args.reverse_convert_out)
