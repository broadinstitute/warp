import hail as hl
import argparse
import os

def write_vcf(inputs):
    OutputBucket = inputs['OutputBucket'] 
    if not OutputBucket.endswith('/'):
        OutputBucket += '/'
    
    OutputFilePath = OutputBucket + inputs['OutputPrefix'] + '.vcf.bgz'
    print('Writing VCF to:')
    print(OutputFilePath)

    #LOAD TABLES AND FIND SUBSET
    print('hail reading matrix table')
    mt = hl.read_matrix_table(inputs['MatrixTable'])
    
    hl.export_vcf(mt, OutputFilePath)
    print('Finished writing VCF')
    
    with open("outpath.txt", "w") as f:
        f.write(OutputFilePath)   

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--MatrixTable", required=True)
    parser.add_argument("--OutputBucket", required=True)
    parser.add_argument("--OutputPrefix", required=True)
    parser.add_argument("--CloudTmpdir", required=True)


    args = parser.parse_args()

    inputs = {
        'MatrixTable': args.MatrixTable,
        'OutputBucket': args.OutputBucket,
        'OutputPrefix': args.OutputPrefix
    }

    hl.init(
        app_name='hail_job',
        master='local[*]',
        tmp_dir=f'{args.CloudTmpdir}',  # Cloud storage recommended here
        spark_conf={
            'spark.executor.instances': '4',
            'spark.executor.cores': '8',
            'spark.executor.memory': '25g',
            'spark.driver.memory': '30g',
            'spark.local.dir': '/cromwell_root',
            'spark.sql.shuffle.partitions': '100',
            'spark.default.parallelism': '100',
            'spark.memory.fraction': '0.8',
            'spark.memory.storageFraction': '0.2',
        },
        default_reference='GRCh38'
    )
    
    print("Spark local directories:", os.getenv("SPARK_LOCAL_DIRS"), flush=True)
    print("Disk usage:", flush=True)
    os.system("df -h")
    
    write_vcf(inputs)
    
