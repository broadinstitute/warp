import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Remove intervals that blew up after liftover')
parser.add_argument('input')
parser.add_argument('output')

args = parser.parse_args()
df = pd.read_csv(args.input, sep='\t', names=("chr", "start", "end", "strand", "name"))
df[['const', 'orig_chr', 'orig_start', 'orig_end']] = df['name'].str.split('-', expand=True)
df = df.astype({'orig_start': 'int32', 'orig_end': 'int32'})
good = df[(df["end"]-df["start"])/(df["orig_end"]-df["orig_start"]) < 10]
out = good[['chr', 'start', 'end', 'strand', 'name']]
out.to_csv(args.output, sep='\t', header=False, index=False)

