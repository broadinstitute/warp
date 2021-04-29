#!/usr/bin/env python3
import argparse
import gzip
import re
import sys
from bisect import bisect_left, bisect_right


def get_feature(line, feature):
    features = re.sub('\"', '', line.strip().split('\t')[8].strip())
    features_dic = {x.split()[0]:x.split()[1] for x in features.split(';') if x} 

    if feature in features_dic:
       return features_dic[feature]
    return None
    

def main():
    """ This script adds intronic features to a GTF file for single-nuclei processing.
        and subsequently use it with featurecounts to show intronic counts
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
      "--input-gtf", "-i", dest="input_gtf", default=None, required=True, help="input GTF"
    )
    parser.add_argument(
      "--output-gtf", "-o", dest="output_gtf", default=None, help="output GTF"
    )
    args = parser.parse_args()

    writeBinary = False


    intron_cands = {}
    transcript_gene = {}

    exons = {}
    exon_ids = {}
    gene_locs = {}
    gene_locations ={}
    # gather the location of each genes and exons; and exon_ids to avoid duplication
    with gzip.open(args.input_gtf, "rt") if args.input_gtf.endswith(".gz"
    ) else open(args.input_gtf, "r") as input_file:
        for line in input_file:
            if not line.startswith("#"):
               fields = [x.strip() for x in line.strip().split('\t')]
               if fields[2] == 'exon':
                  gene_id = get_feature(line.strip(), 'gene_id')
                  exon_id = get_feature(line.strip(), 'exon_id')
                  contig_id = fields[0] 
                  locpair = (int(fields[3]), int(fields[4]))
                  if contig_id not in exons:
                     exons[contig_id] = []
                  if exon_id not in exon_ids:
                      exons[contig_id].append(locpair)
                      exon_ids[exon_id] = True
               elif fields[2] == 'gene':
                  gene_id = get_feature(line.strip(), 'gene_id')
                  contig_id = fields[0] 
                  locpair = (int(fields[3]), int(fields[4]), gene_id)
                  if gene_id!=None:    
                     if contig_id not in gene_locs:
                         gene_locs[contig_id] = []
                     gene_locs[contig_id].append(locpair)

                     gene_locations[gene_id] = locpair

    
    # sorted the gene locs by start
    for contig_id in gene_locs:
        gene_locs[contig_id].sort(key = lambda x: x[0], reverse=False)
        #print(contig_id, len(gene_locs[contig_id]), gene_locs[contig_id][:3], gene_locs[contig_id][-3:])

    # keep sort the exons by start by contig
    for contig_id in exons:
        exons[contig_id].sort(key = lambda x: x[0], reverse=False)

    # compute the intron candidates for each contig
    # where any bp that is not an exon is an candidate intron whithout 
    # worrying about the inclusiveness of that base pair within the range
    # of a gene
    for contig_id in exons:
        intron_cands[contig_id] = []
        last_exon_end = 0
        for exon_coor in exons[contig_id]:
            if exon_coor[0] > last_exon_end:
               pair = (last_exon_end, exon_coor[0])
               intron_cands[contig_id].append(pair)
               if contig_id=='chr7' and not(pair[1] <= gene_locations['ENSMUSG00000041141'][0] or 
                     pair[0] >= gene_locations['ENSMUSG00000041141'][1]):  # chr7, ENSMUSG00000041141
                   pass
                   #print(contig_id, 'intron candid',  pair[0], pair[1], ' gene loc ', gene_locations['ENSMUSG00000041141'])
         


            last_exon_end = max(last_exon_end, exon_coor[1])

        #add the remaining last  
        pair = (last_exon_end, 30000000000)
        intron_cands[contig_id].append(pair)

        if contig_id=='chr7' and not(pair[1] <= gene_locations['ENSMUSG00000041141'][0] or 
                pair[0] >= gene_locations['ENSMUSG00000041141'][1]):  # chr7, ENSMUSG00000041141
           pass
           #print(contig_id, 'intron candid',  pair[0], pair[1], ' gene loc ', gene_locations['ENSMUSG00000041141'])
         
    
    #sys.exit(0)
    #[ (2, 3), (4, 5), (34, 78) ]
    #[ 2, 3, 4, 5 ]
    # global ordered (ascending) arry of intronic start or end points

    introns = {}
    for contig_id in gene_locs:

        introns[contig_id] = []
        intronic_points = []
        for coor in intron_cands[contig_id]:
            intronic_points.append(coor[0])
            intronic_points.append(coor[1])

        for gene_loc in gene_locs[contig_id]:
           i =  bisect_right(intronic_points, gene_loc[0], 0, len(intronic_points))
           j =  bisect_left(intronic_points, gene_loc[1], 0, len(intronic_points))

           #print(gene_loc, intronic_points[i:j + 1], [x%2 for x in range(i, j + 1)])
            
           if i%2 == 1: # it is a start location on i
              intron_start = gene_loc[0]
              intron_end = intronic_points[i]
              #introns[contig_id].append( (intron_start, intron_end, gene_loc[2]) )
       
              if contig_id=='chr7' and not(intron_end <= gene_locations['ENSMUSG00000041141'][0] or 
                  intron_start >= gene_locations['ENSMUSG00000041141'][1]):  # chr7, ENSMUSG00000041141
                  pass
                  #print(contig_id, 'intron',  intron_start, intron_end, ' gene loc ', gene_locations['ENSMUSG00000041141'])
         
           for k in range(i, j, 2):
              introns[contig_id].append( (intronic_points[k], intronic_points[k+1], gene_loc[2]) )
              if contig_id=='chr7' and not(intronic_points[k] <= gene_locations['ENSMUSG00000041141'][0] or 
                  intronic_points[k+1] >= gene_locations['ENSMUSG00000041141'][1]):  # chr7, ENSMUSG00000041141
                  pass
                  #print(contig_id, 'intron', intronic_points[k], intronic_points[k+1], ' gene loc ', gene_locations['ENSMUSG00000041141'])

           if j%2 == 1:
              intron_start = intronic_points[j]
              intron_end = gene_loc[1]
              introns[contig_id].append((intron_start, intron_end, gene_loc[2]))
              if contig_id=='chr7' and not(intron_end <= gene_locations['ENSMUSG00000041141'][0] or 
                  intron_start >= gene_locations['ENSMUSG00000041141'][1]):  # chr7, ENSMUSG00000041141
                  pass
                  #print(contig_id, 'intron',  intron_start, intron_end, ' gene loc ', gene_locations['ENSMUSG00000041141'])
      
           #print(intronic_points[i], "({})".format(i), intronic_points[i+1], gene_loc, 
           #      intronic_points[j], "({})".format(j),  intronic_points[j+1]) 
            
    #sys.exit(0)
    genewise_introns = {}
    for contig_id in introns:
        genewise_introns[contig_id] = {}
        for intron in introns[contig_id]:
            if intron[2] not in genewise_introns[contig_id]:
                genewise_introns[contig_id][intron[2]] = []
            genewise_introns[contig_id][intron[2]].append((intron[0], intron[1]))
              
         #print(contig_id, len(introns[contig_id]), introns[contig_id][:5])
    intron_no = 1
    with gzip.open(args.input_gtf, "rt") if args.input_gtf.endswith(".gz"
    ) else open(args.input_gtf, "r") as input_file:
        for line in input_file:
            if line.startswith("#"):
               print(line.strip())
            else:
               fields = [x.strip() for x in line.strip().split('\t')]
               if fields[2] == 'exon':
                  print(line.strip())
               elif fields[2] == 'gene':
                   print(line.strip())
                   gene_id = get_feature(line.strip(), 'gene_id')
                   contig_id = fields[0]
                   if gene_id in genewise_introns[contig_id]:
                       for intron in genewise_introns[contig_id][gene_id]:
                           mod_fields = fields.copy()
                           mod_fields[2] = 'intron'
                           mod_fields[3] = str(intron[0])
                           mod_fields[4] = str(intron[1])
                           mod_fields[8] = mod_fields[8] + " intron_id \"{}\"".format(str(intron_no))
                           intron_no += 1
                           print('\t'.join(mod_fields))
               else:
                   print(line.strip())

    sys.exit(0)

if __name__ == "__main__":
    main()
