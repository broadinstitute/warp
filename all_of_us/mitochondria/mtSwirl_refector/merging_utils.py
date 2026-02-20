import math
import os
import hail as hl
import logging

from copy import deepcopy
from typing import Dict
from hail.utils.java import info
from generate_mtdna_call_mt.merging_constants import *


def multi_way_union_mts(mts: list, temp_dir: str, chunk_size: int, min_partitions: int, check_from_disk: bool, prefix: str) -> hl.MatrixTable:
    """
    Hierarchically join together MatrixTables in the provided list.

    :param mts: List of MatrixTables to join together
    :param temp_dir: Path to temporary directory for intermediate results
    :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
    :return: Joined MatrixTable
    """
    # Convert the MatrixTables to tables where entries are an array of structs
    if check_from_disk:
        staging = [x for x in mts]
    else:
        staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
    
    stage = 0
    while len(staging) > 1:
        # Calculate the number of jobs to run based on the chunk size
        n_jobs = int(math.ceil(len(staging) / chunk_size))
        info(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
        
        next_stage = []
        
        if check_from_disk:
            all_exists = True
            for idx in range(n_jobs):
                path = os.path.join(temp_dir, f"stage_{stage}_job_{idx}.ht")
                exists = hl.hadoop_is_file(f'{path}/_SUCCESS')
                if not exists:
                    print(path + ' is missing.')
                    if stage == 0:
                        raise ValueError('ERROR: --check-from-disk was enabled but not all stage 0 MTs were found. This is unsupported.')
                    all_exists = False
                    break
            
            if all_exists:
                info(f"Reading stage {stage} from disk...")
                staging.clear()
                for idx in range(n_jobs):
                    staging.append(hl.read_table(os.path.join(temp_dir, f"stage_{stage}_job_{idx}.ht")))
                info(f"Stage {stage} imported from disk.")
                stage += 1
                continue

        for i in range(n_jobs):
            # Grab just the tables for the given job
            to_merge = staging[chunk_size * i : chunk_size * (i + 1)]
            info(
                f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs"
            )

            # Multiway zip join will produce an __entries annotation, which is an array where each element is a struct containing the __entries annotation (array of structs) for that sample
            merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
            if min_partitions > 10:
                merged = merged.repartition(min_partitions).checkpoint(os.path.join(temp_dir, f"{prefix}stage_{stage}_job_{i}_pre.ht"), overwrite=True)
            # Flatten __entries while taking into account different entry lengths at different samples/variants (samples lacking a variant will be NA)
            merged = merged.annotate(
                __entries=hl.flatten(
                    hl.range(hl.len(merged.__entries)).map(
                        # Coalesce will return the first non-missing argument, so if the entry info is not missing, use that info, but if it is missing, create an entries struct with the correct element type for each null entry annotation (such as int32 for DP)
                        lambda i: hl.coalesce(
                            merged.__entries[i].__entries,
                            hl.range(hl.len(merged.__cols[i].__cols)).map(
                                lambda j: hl.null(
                                    merged.__entries.__entries.dtype.element_type.element_type
                                )
                            ),
                        )
                    )
                )
            )

            # Flatten col annotation from array<struct{__cols: array<struct{s: str}>} to array<struct{s: str}>
            merged = merged.annotate_globals(
                __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
            )

            next_stage.append(
                merged.checkpoint(
                    os.path.join(temp_dir, f"{prefix}stage_{stage}_job_{i}.ht"), overwrite=True
                )
            )
        info(f"Completed stage {stage}")
        stage += 1
        staging.clear()
        staging.extend(next_stage)

    # Unlocalize the entries, and unfilter the filtered entries and populate fields with missing values
    return (
        staging[0]
        ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
        .unfilter_entries()
    )


def chunks(items, binsize):
    lst = []
    for item in items:
        lst.append(item)
        if len(lst) == binsize:
            yield lst
            lst = []
    if len(lst) > 0:
        yield lst


def coverage_merging(paths, num_merges, chunk_size, check_from_disk, 
                     temp_dir, n_read_partitions, n_final_partitions, 
                     keep_targets, logger, no_batch_mode=False):
    
    if no_batch_mode:
        pairs_for_coverage = paths.annotate(pairs = (paths.s, paths.coverage)).pairs.collect()
    else:
        pairs_for_coverage = paths.annotate(pairs = (paths.batch, paths.coverage)).pairs.collect()
    
    if num_merges > 1:
        # check_from_disk is not compatible with multiple merges and will not be used
        merged_prefix = f'coverage_merging_final_{str(num_merges)}subsets/'
        this_merged_mt = os.path.join(temp_dir, f"{merged_prefix}final_merged.mt")
        if hl.hadoop_is_file(this_merged_mt + '/_SUCCESS'):
            cov_mt = hl.read_matrix_table(this_merged_mt)
        else:
            subsets = chunks(pairs_for_coverage, len(pairs_for_coverage) // num_merges)
            mt_list_subsets = []
            for subset_number, subset in enumerate(subsets):
                print(f'Importing subset {str(subset_number)}...')
                this_prefix = f'coverage_merging_subset{str(subset_number)}_{str(num_merges)}subsets/'
                this_subset_mt = os.path.join(temp_dir, f"{this_prefix}final_merged.mt")
                if hl.hadoop_is_file(f'{this_subset_mt}/_SUCCESS'):
                    mt_list_subsets.append(hl.read_matrix_table(this_subset_mt))
                    print(f'Subset {str(subset_number)} already processed and imported with {str(mt_list_subsets[len(mt_list_subsets)-1].count_cols())} samples.')
                else:
                    mt_list = []
                    idx = 0
                    for batch, base_level_coverage_metrics in subset:
                        idx+=1
                        mt = hl.import_matrix_table(
                            base_level_coverage_metrics,
                            delimiter="\t",
                            row_fields={"chrom": hl.tstr, "pos": hl.tint, "target": hl.tstr},
                            row_key=["chrom", "pos"],
                            min_partitions=n_read_partitions,
                        )
                        if not keep_targets:
                            mt = mt.drop("target")
                        else:
                            mt = mt.key_rows_by(*["chrom", "pos", "target"])

                        if no_batch_mode:
                            mt = mt.key_cols_by().annotate_cols(col_id = batch)
                            mt = mt.rename({"x": "coverage", 'col_id':'s'}).key_cols_by('s')
                        else:
                            mt = mt.key_cols_by().rename({"x": "coverage", 'col_id':'s'}).key_cols_by('s')
                            mt = mt.annotate_cols(batch = batch)

                        mt_list.append(mt)
                        if idx % 10 == 0:
                            logger.info(f"Imported batch {str(idx)}, subset {str(subset_number)}...")

                    logger.info(f"Joining individual coverage mts for subset {str(subset_number)}...")
                    cov_mt_this = multi_way_union_mts(mt_list, temp_dir, chunk_size, min_partitions=n_read_partitions, check_from_disk=False, prefix=this_prefix)
                    cov_mt_this = cov_mt_this.repartition(n_final_partitions // num_merges).checkpoint(this_subset_mt, overwrite=True)
                    mt_list_subsets.append(cov_mt_this)
            cov_mt = multi_way_union_mts(mt_list_subsets, temp_dir, chunk_size, min_partitions=n_read_partitions, check_from_disk=False, prefix=merged_prefix)
            cov_mt = cov_mt.repartition(n_final_partitions).checkpoint(this_merged_mt, overwrite=True)
    else:
        mt_list = []
        idx = 0
        if check_from_disk:
            logger.info("NOTE: Skipping reading individual coverage MTs since --check-from-disk was enabled.")
            n_append = len(pairs_for_coverage)-1
            pairs_for_coverage = pairs_for_coverage[0:1]
        
        for batch, base_level_coverage_metrics in pairs_for_coverage:
            idx+=1
            mt = hl.import_matrix_table(
                base_level_coverage_metrics,
                delimiter="\t",
                row_fields={"chrom": hl.tstr, "pos": hl.tint, "target": hl.tstr},
                row_key=["chrom", "pos"],
                min_partitions=n_read_partitions,
            )
            if not keep_targets:
                mt = mt.drop("target")
            else:
                mt = mt.key_rows_by(*["chrom", "pos", "target"])
            
            if no_batch_mode:
                mt = mt.key_cols_by().annotate_cols(col_id = batch)
                mt = mt.rename({"x": "coverage", 'col_id':'s'}).key_cols_by('s')
            else:
                mt = mt.key_cols_by().rename({"x": "coverage", 'col_id':'s'}).key_cols_by('s')
                mt = mt.annotate_cols(batch = batch)

            mt_list.append(mt)
            if idx % 10 == 0:
                logger.info(f"Imported batch {str(idx)}...")

        if check_from_disk:
            mt_list.extend([None for x in range(n_append)])

        logger.info("Joining individual coverage mts...")
        cov_mt = multi_way_union_mts(mt_list, temp_dir, chunk_size, min_partitions=n_read_partitions, check_from_disk=check_from_disk, prefix='')

    cov_mt = cov_mt.annotate_rows(locus = hl.locus(cov_mt.chrom, cov_mt.pos, reference_genome="GRCh38"))
    cov_mt = cov_mt.key_rows_by("locus").drop("chrom", "pos")
    return cov_mt


def join_two_mts(mt1, mt2, row_keep, col_keep, temp_dir, partitions):

    mt1 = mt1.select_rows(*row_keep).select_cols(*col_keep)
    mt2 = mt2.select_rows(*row_keep).select_cols(*col_keep)

    if (mt1.row.dtype != mt2.row.dtype) or (mt1.col.dtype != mt2.col.dtype) or (mt1.entry.dtype != mt2.entry.dtype):
        raise ValueError('ERROR: when joining MatrixTables, schemas must be the same.')
    
    mt_append = multi_way_union_mts(mts=[mt1,mt2], temp_dir=temp_dir, chunk_size=2, min_partitions=partitions, check_from_disk=False, prefix='appending_')
    return mt_append


def append_coverage_to_old(cov_mt, old_mt_path, col_keep, n_final_partitions, temp_dir):

    this_merged_mt = os.path.join(temp_dir, 'coverage_tmp_appended_to_old_dataset_final.mt')
    cov_mt = cov_mt.checkpoint(os.path.join(temp_dir, 'coverage_mt_new_keyed_pre_merge_with_old.mt'), overwrite=True)

    if hl.hadoop_is_file(this_merged_mt + '/_SUCCESS'):
        cov_mt = hl.read_matrix_table(this_merged_mt)
    else:
        old_mt = hl.read_matrix_table(old_mt_path)
        print(f'Second database imported with {str(old_mt.count()[1])} samples.')
        cov_mt = join_two_mts(mt1 = old_mt, mt2 = cov_mt, row_keep = [], col_keep = col_keep, temp_dir=temp_dir, partitions=n_final_partitions)
        cov_mt = cov_mt.repartition(n_final_partitions).checkpoint(this_merged_mt, overwrite=True) 

    return cov_mt


def add_coverage_annotations(cov_mt):
    n_samples = cov_mt.count_cols()
    cov_mt = cov_mt.annotate_rows(
        mean=hl.float(hl.agg.mean(cov_mt.coverage)),
        median=hl.median(hl.agg.collect(cov_mt.coverage)),
        over_100=hl.float((hl.agg.count_where(cov_mt.coverage > 100) / n_samples)),
        over_1000=hl.float((hl.agg.count_where(cov_mt.coverage > 1000) / n_samples)),
    )
    return cov_mt


def append_vcf_to_old(combined_mt, old_mt_path, col_keep, n_final_partitions, temp_dir):
    this_merged_mt = os.path.join(temp_dir, 'variants_tmp_appended_to_old_dataset_final.mt')

    if hl.hadoop_is_file(this_merged_mt + '/_SUCCESS'):
        combined_mt = hl.read_matrix_table(this_merged_mt)
    else:
        old_mt = hl.read_matrix_table(old_mt_path)
        print(f'Second database imported with {str(old_mt.count()[1])} samples and {str(old_mt.count()[0])} variants.')
        combined_mt = join_two_mts(mt1 = old_mt, mt2 = combined_mt, row_keep = [], col_keep = col_keep, temp_dir=temp_dir, partitions=n_final_partitions)
        combined_mt = combined_mt.repartition(n_final_partitions).checkpoint(this_merged_mt, overwrite=True) 

    return combined_mt


def vcf_merging_and_processing(
    vcf_paths,
    coverage_mt_path,
    include_extra_v2_fields,
    single_sample,
    old_mt_path,
    artifact_prone_sites_path,
    artifact_prone_sites_reference,
    minimum_homref_coverage,
    logger,
    chunk_size,
    num_merges,
    n_final_partitions,
    output_bucket,
    temp_dir,
    overwrite,
    *,
    coverage_h5_path: str | None = None,
    homref_position_block_size: int = 1024,
):
    """v1/v2 compatible wrapper.

    If `coverage_h5_path` is provided, use covdb-based hom-ref imputation.
    Otherwise fall back to v1 behavior using `coverage_mt_path`.
    """
    output_path_mt = os.path.join(output_bucket, "raw_combined.mt")
    output_path_mt_2 = os.path.join(output_bucket, "raw_combined_2.mt")

    if hl.hadoop_exists(f"{output_path_mt}/_SUCCESS") and not overwrite:
        logger.info(f"Reading merged VCF mt from {output_path_mt}...")
        combined_mt = hl.read_matrix_table(output_path_mt, _n_partitions=n_final_partitions)
        meta = get_vcf_metadata(include_extra_v2_fields)
    else:
        logger.info("Combining VCFs...")
        combined_mt, meta = vcf_merging(
            vcf_paths=vcf_paths,
            temp_dir=temp_dir,
            logger=logger,
            chunk_size=chunk_size,
            include_extra_v2_fields=include_extra_v2_fields,
            num_merges=num_merges,
            single_sample=single_sample,
            n_final_partitions=n_final_partitions,
        )
        combined_mt = combined_mt.naive_coalesce(n_final_partitions).checkpoint(
            output_path_mt, overwrite=overwrite
        )
        combined_mt = hl.read_matrix_table(output_path_mt, _n_partitions=n_final_partitions)

    if hl.hadoop_exists(f"{output_path_mt_2}/_SUCCESS") and not overwrite:
        logger.info(f"Reading merged VCF mt from {output_path_mt_2}...")
        combined_mt = hl.read_matrix_table(output_path_mt_2, _n_partitions=n_final_partitions)
    else:
        logger.info("Removing select sample-level filters...")
        combined_mt = remove_genotype_filters(combined_mt)

        logger.info("Determining homoplasmic reference sites...")
        if coverage_h5_path is not None:
            combined_mt = determine_hom_refs_from_covdb(
                combined_mt,
                coverage_h5_path=coverage_h5_path,
                minimum_homref_coverage=minimum_homref_coverage,
                position_block_size=homref_position_block_size,
                logger=logger,
            )
        else:
            combined_mt = determine_hom_refs(
                combined_mt, coverage_mt_path, minimum_homref_coverage
            )
        combined_mt = combined_mt.checkpoint(output_path_mt_2, overwrite=overwrite)

    if old_mt_path is not None:
        logger.info("Appending new VCF to old VCF database...")
        col_keep = [] if single_sample else ["batch"]
        combined_mt = append_vcf_to_old(
            combined_mt, old_mt_path, col_keep, n_final_partitions, temp_dir
        )

    logger.info("Applying artifact_prone_site filter...")
    combined_mt = apply_mito_artifact_filter(
        combined_mt, artifact_prone_sites_path, artifact_prone_sites_reference
    )

    return combined_mt, meta


def get_vcf_metadata(include_extra_v2_fields: bool):
    meta = deepcopy(META_DICT_BASE)
    if include_extra_v2_fields:
        meta['format'].update(META_DICT_V2_FMT)

    return meta


def vcf_merging(vcf_paths: Dict[str, str], temp_dir: str, logger, n_final_partitions, 
                chunk_size: int = 100, include_extra_v2_fields: bool = False, num_merges: int = 1,
                single_sample: bool = False) -> hl.MatrixTable:
    """
    Reformat and join individual mitochondrial VCFs into one MatrixTable.

    :param vcf_paths: Dictionary of samples to combine (sample as key, path to VCF as value)
    :param temp_dir: Path to temporary directory for intermediate results
    :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
    :param include_extra_v2_fields: Includes extra fields important for analysis of v2.1 source MTs
    :return: Joined MatrixTable of samples given in vcf_paths dictionary
    """
    # Update VCF metadata
    meta = get_vcf_metadata(include_extra_v2_fields)

    list_paths = list(vcf_paths.items())
    list_paths.sort(key=lambda y: y[0])
    if num_merges == 1:
        vcf_path_list = [list_paths]
    else:
        vcf_path_list = chunks(list_paths, len(list_paths) // num_merges)
    mt_list_subsets = []
    for subset_number, subset in enumerate(vcf_path_list):
        print(f'Importing subset {str(subset_number)}...')
        this_prefix = f'variant_merging_subset{str(subset_number)}_{str(num_merges)}subsets/'
        this_subset_mt = os.path.join(temp_dir, f"{this_prefix}final_merged.mt")
        if hl.hadoop_is_file(f'{this_subset_mt}/_SUCCESS'):
            mt_list_subsets.append(hl.read_matrix_table(this_subset_mt))
            print(f'Subset {str(subset_number)} already processed and imported with {str(mt_list_subsets[len(mt_list_subsets)-1].count_cols())} samples.')
        else:
            mt_list = []
            idx = 0
            for batch, vcf_path in subset:
                idx+=1
                try:
                    mt = hl.import_vcf(vcf_path, reference_genome="GRCh38", array_elements_required=False)
                except Exception as e:
                    raise ValueError(
                        f"vcf path {vcf_path} does not exist for sample {batch}"
                    ) from e

                # Because the vcfs are split, there is only one AF value, although misinterpreted as an array because Number=A in VCF header
                # Second value of MMQ is the value of the mapping quality for the alternate allele
                # Add FT annotation for sample genotype filters (pull these from filters annotations of the single-sample VCFs)
                if include_extra_v2_fields:
                    if 'GT' in mt.entry:
                        mt = mt.drop('GT')
                    
                if single_sample:
                    if include_extra_v2_fields:
                        # process keys that are already entries
                        entry_fields = {k:v for k, v in V2_FIELD_KEY.items() if k not in V2_INFO_TO_FORMAT}
                        for x, item_type in entry_fields.items():
                            if x not in mt.entry:
                                mt = mt.annotate_entries(**{x: hl.missing(item_type)})
                        mt = mt.select_entries("DP", "AD", *list(entry_fields.keys()), HL=mt.AF[0])

                        # now process info fields into entries
                        info_fields = {k:v for k, v in V2_FIELD_KEY.items() if k in V2_INFO_TO_FORMAT}
                        for x, item_type in info_fields.items():
                            if x not in mt.info:
                                mt = mt.annotate_entries(**{x: hl.missing(item_type)})
                            else:
                                if x in V2_INFO_TO_FORMAT_REQUIREINDEX:
                                    mt = mt.annotate_entries(**{x: mt.info[x][0]})
                                else:
                                    mt = mt.annotate_entries(**{x: mt.info[x]})
                                if mt.info[x].dtype == hl.dtype('tbool'):
                                    # flag is not supported in FORMAT
                                    mt = mt.annotate_entries(**{x: hl.if_else(mt.info[x], 1, 0)})
                    else:
                        mt = mt.select_entries("DP", HL=mt.AF[0])

                    mt = mt.annotate_entries(
                        MQ=hl.float(mt.info["MMQ"][1]),
                        TLOD=mt.info["TLOD"][0],
                        FT=hl.if_else(hl.len(mt.filters) == 0, {"PASS"}, mt.filters)
                    )
                else:
                    if include_extra_v2_fields:
                        for x, item_type in V2_FIELD_KEY.items():
                            if x not in mt.entry:
                                mt = mt.annotate_entries(**{x: hl.missing(item_type)})
                        
                # enforce entry ordering
                if include_extra_v2_fields:
                    mt = mt.select_entries("DP", "AD", *list(V2_FIELD_KEY.keys()), "HL", "MQ", "TLOD", "FT")
                else:
                    mt = mt.select_entries("DP", "HL", "MQ", "TLOD", "FT")
                
                # Use GRCh37 reference as most external resources added in downstream scripts use GRCh37 contig names
                # (although note that the actual sequences of the mitochondria in both GRCh37 and GRCh38 are the same)
                mt = mt.annotate_entries(FT = hl.set(mt.FT))
                mt = mt.key_rows_by(
                    locus=hl.locus("MT", mt.locus.position, reference_genome="GRCh37"),
                    alleles=mt.alleles,
                )
                if single_sample:
                    mt = mt.key_cols_by(s=batch)
                else:
                    mt = mt.annotate_cols(batch=batch).key_cols_by('s')
                mt = mt.select_rows()
                mt_list.append(mt)
                if idx % 20 == 0:
                    if single_sample:
                        logger.info(f"Imported sample {str(idx)}...")
                    else:    
                        logger.info(f"Imported batch {str(idx)}...")

            n_indiv_repart = n_final_partitions // num_merges
            combined_mt_this = multi_way_union_mts(mt_list, temp_dir, chunk_size, min_partitions=n_indiv_repart // 2, check_from_disk=False, prefix=this_prefix)
            combined_mt_this = combined_mt_this.repartition(n_indiv_repart).checkpoint(this_subset_mt, overwrite=True)
            mt_list_subsets.append(combined_mt_this)
    
    if num_merges == 1:
        combined_mt = mt_list_subsets[0]
    else:
        merged_prefix = f'variant_merging_final_{str(num_merges)}subsets/'
        combined_mt = multi_way_union_mts(mt_list_subsets, temp_dir, chunk_size, min_partitions=n_indiv_repart // 2, check_from_disk=False, prefix=merged_prefix)

    return combined_mt, meta


def collect_vcf_paths(participant_data: str, vcf_col_name: str, participants_to_subset: str = None,
                      single_sample: bool = False) -> Dict[str, str]:
    """
    Create dictionary of VCF paths for only the samples specified in participants_to_subset.

    .. note::
        Participant data should be a tab-delimited file with at minimum columns for:
        - 'entity:participant_id': sample name with prohibited characters replaced with underscores
        - 's': sample name
        - path to the Mutect2 VCF output, where name of this column is supplied to the `vcf_col_name` parameter

    :param participant_data: Participant data (a ht)
    :param vcf_col_name: Name of column that contains VCF output
    :param participants_to_subset: Path to file of participant_ids to which the data should be subset
    :return: Dictionary with sample name as key and path to VCF as value
    """
    vcf_paths = {}
    
    # Load in data
    if os.path.splitext(participant_data)[1] == '.ht':
        participant_ht = hl.read_table(participant_data)
    else:
        participant_ht = hl.import_table(participant_data)
    
    # Remove participants that don't have VCF output
    participant_ht.filter(participant_ht[vcf_col_name] != "")

    # Subset participants if specified
    if participants_to_subset:
        participants_of_interest = hl.import_table(
            participants_to_subset
        ).participant.collect()
        participant_ht = participant_ht.filter(
            hl.literal(participants_of_interest).contains(
                participant_ht["entity:participant_id"]
            )
        )

    # Add the vcf path to a dictionary with batch name as key
    df = participant_ht.to_pandas()

    for _, row in df.iterrows():
        if single_sample:
            vcf_paths[row["s"]] = row[vcf_col_name]
        else:
            vcf_paths[row["batch"]] = row[vcf_col_name]

    return vcf_paths


def remove_genotype_filters(mt: hl.MatrixTable,
    filters_to_remove: set = {
        "possible_numt",
        "mt_many_low_hets",
        "FAIL",
        "blacklisted_site",
    },
) -> hl.MatrixTable:
    """
    Remove unneeded sample-level genotype filters (in FT field of the VCF) specified by the filters_to_remove parameter.

    By default, remove the 'possible_numt', 'mt_many_low_hets', and 'FAIL' filters because these filters were found to have low performance.
    Also remove the 'blacklisted_site' filter because this filter did not always behave as expected in early GATK versions. This filter can be reimplemented with the apply_mito_artifact_filter function.

    NOTE: This function does not modify any row or column-level features. Thus it is safe
          to run this before we join the VCF mt with another one.

    :param mt: MatrixTable containing genotype filters in the FT field of the VCF that should be removed
    :param filters_to_remove: List of genptype filters (in FT field of VCF) that should be removed from the entries
    :return: MatrixTable with specific genotype filters (in FT field of VCF) removed
    """
    mt = mt.annotate_entries(FT=mt.FT.difference(filters_to_remove))

    # If no filters exist after removing those specified above, set the FT field to PASS
    mt = mt.annotate_entries(FT=hl.if_else(hl.len(mt.FT) == 0, {"PASS"}, mt.FT))

    return mt


def determine_hom_refs(mt: hl.MatrixTable, coverage_mt_path: str, minimum_homref_coverage: int = 100) -> hl.MatrixTable:
    """
    Use coverage to distinguish between homref and missing sites.
    NOTE: This function does not modify any row or column-level features. Thus it is safe
          to run this before we join the VCF mt with another one.

    :param mt: MatrixTable from initial multi-sample merging, without homref sites determined
    :param coverage_mt_path: MatrixTable of sample level coverage at each position (per-sample and per-base; can be generated by running annotate_coverage.py)
    :param minimum_homref_coverage: Minimum depth of coverage required to call a genotype homoplasmic reference rather than missing
    :return: MatrixTable with missing genotypes converted to homref depending on coverage
    """
    # Convert coverage to build GRCh37 to match contig names
    # Note: the mitochondrial reference genome is the same for GRCh38 and GRCh37
    coverages = hl.read_matrix_table(coverage_mt_path)
    coverages = coverages.key_rows_by(
        locus=hl.locus("MT", coverages.locus.position, reference_genome="GRCh37")
    )

    mt = mt.annotate_entries(
        DP=hl.if_else(hl.is_missing(mt.HL), coverages[mt.locus, mt.s].coverage, mt.DP)
    )

    hom_ref_expr = hl.is_missing(mt.HL) & (mt.DP > minimum_homref_coverage)

    mt = mt.annotate_entries(
        HL=hl.if_else(hom_ref_expr, 0.0, mt.HL),
        FT=hl.if_else(hom_ref_expr, {"PASS"}, mt.FT),
        DP=hl.if_else(
            hl.is_missing(mt.HL) & (mt.DP <= minimum_homref_coverage),
            hl.null(hl.tint32),
            mt.DP,
        ),
    )

    return mt


def determine_hom_refs_from_covdb(
    mt: hl.MatrixTable,
    *,
    coverage_h5_path: str,
    minimum_homref_coverage: int = 100,
    position_block_size: int = 1024,
    checkpoint_interval_blocks: int = 1,
    logger: logging.Logger | None = None,
) -> hl.MatrixTable:
    """v2 version of `determine_hom_refs` sourcing DP from `coverage.h5`.

    This preserves v1 semantics exactly while avoiding a full (pos × sample) materialization.

    Strategy
    --------
    * Map each MT sample to its covdb row index.
    * Chunk MT positions into position blocks.
    * For each block, read coverage from HDF5 once and broadcast only that block.
    * Update entries for rows in that block, leaving other rows untouched.
    * Periodically checkpoint to avoid constructing an enormous IR from all blocks.
    """
    from generate_mtdna_call_mt.covdb_utils import open_covdb_index, read_covdb_block
    import numpy as np

    log = logger if logger is not None else logging.getLogger(__name__)
    log.info(
        "determine_hom_refs_from_covdb: loading covdb index from %s", coverage_h5_path
    )
    idx = open_covdb_index(coverage_h5_path)

    if "s" not in mt.col:
        raise ValueError("Expected MT to have column key 's'")

    sample_to_index = hl.literal(idx.sample_to_index)
    mt = mt.annotate_cols(__covdb_sample_index=sample_to_index.get(mt.s))
    mt = mt.add_col_index(name="__col_idx")

    ht_cols = mt.cols()
    col_rows = ht_cols.select(
        cov_i=ht_cols.__covdb_sample_index,
        col_i=ht_cols.__col_idx,
    ).collect()
    missing_samples = [r.s for r in col_rows if r.cov_i is None]
    if missing_samples:
        raise ValueError(
            f"{len(missing_samples)} MT samples missing from coverage.h5 /sample_id (first 5: {missing_samples[:5]})"
        )

    n_cols = len(col_rows)
    covdb_row_by_col = np.empty(n_cols, dtype=np.int64)
    for r in col_rows:
        covdb_row_by_col[int(r.col_i)] = int(r.cov_i)

    ht_rows = mt.rows()
    pos_ht = ht_rows.select(pos=hl.int32(ht_rows.locus.position)).distinct().order_by("pos")
    positions = [int(r.pos) for r in pos_ht.collect()]
    if not positions:
        return mt

    pos_to_covdb_col = {p: idx.pos_to_index.get(p) for p in positions}
    missing_pos = [p for p, ci in pos_to_covdb_col.items() if ci is None]
    if missing_pos:
        raise ValueError(
            f"{len(missing_pos)} MT positions missing from coverage.h5 /pos (first 5: {missing_pos[:5]})"
        )

    n_blocks = int(math.ceil(len(positions) / position_block_size))
    log.info(
        "determine_hom_refs_from_covdb: processing %d position blocks (block_size=%d)",
        n_blocks,
        position_block_size,
    )

    pos_to_block = {}
    block_positions: list[list[int]] = []
    for block_id, start in enumerate(range(0, len(positions), position_block_size)):
        block = positions[start : start + position_block_size]
        block_positions.append(block)
        for p in block:
            pos_to_block[int(p)] = int(block_id)

    pos_to_block_hl = hl.literal({int(np.int32(k)): int(np.int32(v)) for k, v in pos_to_block.items()})
    mt = mt.annotate_rows(__pos=hl.int32(mt.locus.position))
    mt = mt.annotate_rows(__block=hl.int32(pos_to_block_hl.get(mt.__pos)))
    mt = mt.repartition(n_blocks, shuffle=True)

    mt = mt.annotate_entries(__cov=hl.missing(hl.tint32))

    for block_id, block in enumerate(block_positions):
        block_id_i32 = int(np.int32(block_id))
        covdb_cols_block = np.array([pos_to_covdb_col[p] for p in block], dtype=np.int64)
        log.info(
            "covdb read block %d/%d: positions[%d..%d] (n_pos=%d)",
            block_id + 1,
            n_blocks,
            block[0],
            block[-1],
            len(block),
        )

        cov_block = read_covdb_block(
            h5_path=coverage_h5_path,
            sample_indices=covdb_row_by_col,
            pos_indices=covdb_cols_block,
        ).astype(np.int32)
        cov_block_hl = hl.literal(cov_block.T.tolist())

        pos_to_offset_hl = hl.literal(
            {int(np.int32(p)): int(np.int32(i)) for i, p in enumerate(block)}
        )
        offset_expr = pos_to_offset_hl.get(mt.__pos)
        cov_expr = hl.if_else(
            hl.is_defined(offset_expr),
            cov_block_hl[hl.int32(offset_expr)][hl.int32(mt.__col_idx)],
            hl.missing(hl.tint32),
        )

        mt = mt.annotate_entries(
            __cov=hl.if_else(mt.__block == block_id_i32, cov_expr, mt.__cov)
        )

        if (
            checkpoint_interval_blocks > 0
            and (block_id + 1) % checkpoint_interval_blocks == 0
            and (block_id + 1) < n_blocks
        ):
            tmp_path = hl.utils.new_temp_file("covdb_block", extension="mt")
            log.info(
                "Checkpointing covdb blocks at %d/%d -> %s",
                block_id + 1,
                n_blocks,
                tmp_path,
            )
            mt = mt.checkpoint(tmp_path, overwrite=True)

    mt = mt.annotate_entries(DP=hl.if_else(hl.is_missing(mt.HL), mt.__cov, mt.DP))

    hom_ref_expr = hl.is_missing(mt.HL) & (mt.DP > minimum_homref_coverage)
    mt = mt.annotate_entries(
        HL=hl.if_else(hom_ref_expr, 0.0, mt.HL),
        FT=hl.if_else(hom_ref_expr, {"PASS"}, mt.FT),
        DP=hl.if_else(
            hl.is_missing(mt.HL) & (mt.DP <= minimum_homref_coverage),
            hl.null(hl.tint32),
            mt.DP,
        ),
    )

    mt = mt.drop("__cov")
    mt = mt.drop("__covdb_sample_index", "__col_idx", "__pos", "__block")
    return mt


def apply_mito_artifact_filter(mt: hl.MatrixTable, artifact_prone_sites_path: str, artifact_prone_sites_reference: str) -> hl.MatrixTable:
    """
    Add in artifact_prone_site filter.

    NOTE: these features DO modify row-level features. We will run this AFTER any
          any potential merging.

    :param mt: MatrixTable to be annotated with artifact_prone_sites filter
    :param artifact_prone_sites_path: Path to BED file of artifact_prone_sites to flag in the filters column
    :return: MatrixTable with artifact_prone_sites filter
    """
    # Apply "artifact_prone_site" filter to any SNP or deletion that spans a known problematic site.
    #
    # Robust contig normalization:
    # The public mtDNA blacklist BED we use is usually hg38/chrM, but callers sometimes pass
    # artifact_prone_sites_reference="default" (or other non-GRCh38 strings). In that case Hail will
    # happily import intervals on contig "chrM" which then blow up later when we build GRCh37/MT
    # intervals (exact failure: "Contig 'chrM' is not in reference genome 'GRCh37'").
    #
    # For mtDNA, GRCh37 vs GRCh38 positions are identical; only contig naming differs. So we
    # ALWAYS normalize any BED that contains contig "chrM" onto GRCh37/MT by rebuilding intervals
    # using start/end positions.
    # IMPORTANT: Don't import this BED against GRCh37/GRCh38 reference genomes up-front.
    # Some references (including "default" in our WDL wiring) effectively resolve to GRCh37,
    # and then importing chrM intervals immediately throws before we get a chance to normalize.
    #
    # Instead, import as a plain table (no reference), detect chrM via string fields, and only
    # then build GRCh37/MT intervals.
    bed_raw = hl.import_table(
        artifact_prone_sites_path,
        no_header=True,
        delimiter="\t",
        types={"f0": hl.tstr, "f1": hl.tint32, "f2": hl.tint32},
    )
    bed_raw = bed_raw.annotate(
        contig=bed_raw.f0,
        start=bed_raw.f1,
        end=bed_raw.f2,
    )

    # Expand BED intervals into a literal set of positions (chrM -> MT normalization by design).
    bed_positions_ht = (
        bed_raw
    .annotate(pos=hl.range(bed_raw.start + 1, bed_raw.end + 1))
    .explode("pos")
    .select("pos")
    .key_by("pos")
    .distinct()
    )
    bed_positions = set(bed_positions_ht.pos.collect())
    bed_pos_set = hl.literal(bed_positions)
    print("Artifact prone sites imported.")

    start_pos = hl.int32(mt.locus.position)
    end_pos = hl.int32(mt.locus.position + hl.len(mt.alleles[0]) - 1)
    span_overlaps = hl.any(
        lambda p: bed_pos_set.contains(p),
        hl.range(start_pos, end_pos + 1),
    )

    mt = mt.annotate_rows(
        filters=hl.if_else(
            bed_pos_set.contains(start_pos)
            | bed_pos_set.contains(end_pos)
            | span_overlaps,
            {"artifact_prone_site"},
            {"PASS"},
        )
    )

    return mt
