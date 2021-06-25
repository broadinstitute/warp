#!/usr/bin/env python
import argparse
import re
import gzip



def main():
    """
        This script take the intronic and exonic counts produced by featureCounts and the reference gtf file to produce a
        four column csv file as follows
        gene_id,gene_name,intron_counts,exon_counts
        ......
        .....
    """
    parser = argparse.ArgumentParser(description="Set the gene_name within a gtf to be equivalent to the values within gene_id.")
    parser.add_argument('--in-gtf-file', dest='in_gtf', help='input gtf file')
    parser.add_argument('--intron-counts', dest='intron_counts', 
                        help='The count file produced by featureCount for the alignments to the intronic regions')
    parser.add_argument('--exon-counts', dest='exon_counts', 
                        help='The count file produced by featureCount for the alignments to  exonic alignments')
    parser.add_argument('--output', '-o',  dest='output', help='output file')
    args = parser.parse_args()

    gene_id_to_gene_name = get_gene_id_to_gene_name(args.in_gtf)

    gene_count_introns = read_gene_id_count(args.intron_counts)
    gene_count_exons = read_gene_id_count(args.exon_counts)

    gene_ids = list(gene_id_to_gene_name.keys())

    with open(args.output, 'w') as fout:
        fout.write("{}\n".format('\t'.join(["gene_id", "gene_name", "introns",
                                            "intron_FPKM",  "exons", "exon_FPKM", 
                                            "intron_length", "exon_length"])))
        
        for gene_id in gene_ids:
            # intron count else 0
            if gene_id in gene_count_introns: 
                length_of_feature = gene_count_introns[gene_id][0]
                intron_count = gene_count_introns[gene_id][1]
                try:
                   intron_fpkm = intron_count/length_of_feature
                except:
                   intron_fpkm = 0
            else:
                intron_count = 0 
                intron_fpkm = 0
                length_of_feature = 0
            intron_length = length_of_feature
    
            # exon count else 0
            if gene_id in gene_count_exons: 
                length_of_feature = gene_count_exons[gene_id][0]
                exon_count = gene_count_exons[gene_id][1]
                try:
                   exon_fpkm = exon_count/length_of_feature
                except:
                   exon_fpkm = 0
            else:
                exon_count = 0 
                exon_fpkm = 0
                length_of_feature = 0
            exon_length = length_of_feature
                 
            fout.write("{}\t{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{}\n".format(gene_id, gene_id_to_gene_name[gene_id], 
                                                                 intron_count, intron_fpkm, 
                                                                 exon_count, exon_fpkm, intron_length, exon_length))

def read_gene_id_count(count_file):
    gene_id_to_count = {}
    with open(count_file, 'r') as fpin:
        for _line in fpin:
            fields = [x.strip() for x in _line.strip().split('\t')]
            if len(fields)==7:
               try:
                  gene_id_to_count[fields[0]] = [int(fields[5]), int(fields[6])]
               except:
                  pass

    return gene_id_to_count 


def get_gene_id_to_gene_name(in_gtf):
    gene_id_to_gene_name = {}
    with gzip.open(in_gtf, 'rt')  if in_gtf.endswith('.gz') else open(in_gtf, 'r') as fpin:
        for _line in fpin:
            line = _line.strip()
            gene_id_search = re.search(r'gene_id ([^;]*);', line)
            gene_name_search = re.search(r'gene_name ([^;]*);', line)
            if gene_id_search and gene_name_search:
                gene_id= re.sub("\"", '', gene_id_search.group(1))
                gene_name= re.sub("\"", '', gene_name_search.group(1))

                if gene_id not in gene_id_to_gene_name:
                   gene_id_to_gene_name[gene_id] = gene_name 

    return gene_id_to_gene_name 

if __name__=='__main__':
    main()
