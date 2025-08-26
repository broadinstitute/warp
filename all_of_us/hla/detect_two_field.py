import sys

# call this script as python [script_name].py hlahd.txt
# this counts the number of A, B, and C alleles called with less than three fields
# and prints this count to stdout

hlahd_file = sys.argv[1]
less_than_three_field_count = 0
with open(hlahd_file, 'r') as file:
    for line in file:
        tokens = line.strip().split('\t')
        gene = tokens[0]

        if gene in ['A', 'B', 'C']:
            for allele in tokens[1:]:
                if allele != '-':   # dash denotes only one allele due to homozygosity
                    num_fields = len(allele.split(':'))
                    if num_fields < 3:
                        less_than_three_field_count += 1

print(less_than_three_field_count)
