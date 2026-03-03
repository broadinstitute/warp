import hail as hl
import argparse

def main(args):
    # Initialize Hail with hard-coded configuration
    hl.init(
    app_name='hail_job',
    master='local[*]',
    tmp_dir=f'{args.CloudTmpdir}',  # Cloud storage recommended here
    spark_conf={
        'spark.local.dir': '/cromwell_root',  # Local SSD for Spark shuffle/spill
        'spark.executor.instances': '4',
        'spark.executor.cores': '16',
        'spark.executor.memory': '64g',
        'spark.driver.memory': '64g',
        'spark.sql.shuffle.partitions': '100',
        'spark.default.parallelism': '100',
        'spark.memory.fraction': '0.8',
        'spark.memory.storageFraction': '0.2'
    },
    default_reference='GRCh38'
    )

    # Load matrix table and samples table
    mt = hl.read_matrix_table(args.MatrixTable)
    samples_ht = hl.import_table(args.SampleList, key='research_id')

    # Filter matrix table to samples in samples_ht
    mt_filtered = mt.filter_cols(hl.is_defined(samples_ht[mt.s]))
    # Filter for biallelic
    mt_filtered = mt_filtered.filter_rows(hl.len(mt_filtered.alleles) == 2)
    # remove empty GT
    mt_filtered = mt_filtered.filter_rows(hl.agg.any(hl.is_defined(mt_filtered.GT)))
    # remove missing AC
    mt_filtered = mt_filtered.filter_rows(~hl.is_missing(mt_filtered.info.AC))
    # filter out bad quals
    mt_filtered = mt_filtered.filter_rows(hl.is_missing(mt_filtered.filters))

    # add variant qc stats
    mt_filtered = hl.variant_qc(mt_filtered)
    # save off total qc
    mt_filtered = mt_filtered.annotate_rows(
        total = mt_filtered.info.annotate(
            ALL_p_value_hwe = mt_filtered.variant_qc.p_value_hwe,
            ALL_p_value_excess_het = mt_filtered.variant_qc.p_value_excess_het
        )
    )
    # recalculate AC/AF/AN on new filtered set
    mt_filtered = mt_filtered.annotate_rows( info = hl.agg.call_stats(mt_filtered.GT, mt_filtered.alleles) )

    # Allele number percentage cutoff
    mt_filtered = mt_filtered.filter_rows(mt_filtered.info.AN >= int(args.AlleleNumberPercentage)/100 * mt_filtered.count_cols() * 2)
    
    # filter by allele count
    mt_filtered = mt_filtered.filter_rows(hl.min(mt_filtered.info.AC) >= int(args.AlleleCount)) 
    
    # save to info field to export to vcf
    mt_filtered = mt_filtered.annotate_rows(
            info = mt_filtered.info.annotate(
                ALL_p_value_hwe = mt_filtered.total.ALL_p_value_hwe,
                ALL_p_value_excess_het = mt_filtered.total.ALL_p_value_excess_het,
                AF = hl.min(mt_filtered.info.AF),
                AC = hl.min(mt_filtered.info.AC),
                AN = mt_filtered.info.AN,
            )
        )
    # get rid of unneeded fields for matrix table save
    mt_filtered = mt_filtered.drop("variant_qc","total")
    
    # Write filtered matrix table to output checkpoint
    mt_filtered.write(f'{args.OutputBucket}/{args.OutputPrefix}.mt', overwrite=True)

    hl.stop()

    with open('outpath.txt', 'w') as file:
        file.write(f'{args.OutputBucket}/{args.OutputPrefix}.mt')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and write Hail MatrixTable with hard-coded Hail configuration.")
    parser.add_argument("--MatrixTable", required=True, help="Path to input MatrixTable.")
    parser.add_argument("--SampleList", required=True, help="Path to samples TSV file.")
    parser.add_argument("--AlleleCount", required=True, help="Allele count threshold.")
    parser.add_argument("--AlleleNumberPercentage", required=True, help="Allele number percentage cutoff.")
    parser.add_argument("--OutputBucket", required=True, help="Path to output checkpoint MatrixTable.")
    parser.add_argument("--OutputPrefix", required=True, help="Output prefix.")
    parser.add_argument("--CloudTmpdir", required=True, help="Temporary directory for spark/hail to work with.")

    args = parser.parse_args()
    main(args)
