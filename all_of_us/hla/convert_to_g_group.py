import sys

# call this script as python [script_name].py hlahd.txt hla_nom_g_groups.txt
# outputs to stdout

hlahd_file = sys.argv[1]
g_group_file = sys.argv[2]
allele_to_group_map = {}

# process the g group file (from https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_g.txt)
# Here are some examples of lines in this file and what they mean

# Example 1:
# A*;11:03:01:01/11:03:01:02/11:175/11:348;11:03:01G
# This means that alleles A*11:03:01:01, A*11:03:01:02, A*11:175, and A*11:348 all
# correspond to 3-field G group allele A*11:03:01 (the token after the last semicolon)

# Example 2:
# A*;11:02:02;
# This is a singleton group containing only one allele.  The allele and group are both A*11:02:02

with open(g_group_file, 'r') as file:
    for line in file:
        # split by semicolons.  The first token is eg A*, the HLA gene.  The second token is a list of alleles,
        # separated by slashes eg 11:03:01:01/11:03:01:02.  The third and last token is the group eg 11:03:01G.
        tokens = line.strip().split(';')

        if tokens[0].startswith('#'):
            # header line, ignore
            continue

        assert len(tokens) == 3
        gene = tokens[0]	# this will be something like A*, for example
        alleles = tokens[1].split('/')
        group = gene + (tokens[-1] if tokens[-1] != '' else alleles[-1])

        # 'G' is considered default group type
        if group.endswith('G'):
            group = group[:-1]

        # three-field alleles are assigned to the group.  Four field alleles are truncated to three fields.
        # two-field alleles are always singletons, whence the allele IS the group, as above.
        for allele in alleles:
            allele_to_group_map[gene + allele] = group
            if allele.count(':') > 2:
                three_field_allele = (':').join(allele.split(':')[:3])
                allele_to_group_map[gene + three_field_allele] = group

        # this handles the singleton two-field allele case
        allele_to_group_map[group] = group

with open(hlahd_file, 'r') as file:
    for line in file:
        tokens = line.strip().split('\t')
        gene = tokens[0]
        new_tokens = [gene]
        for token in tokens[1:]:
            # replace every HLA allele in the file with its group
            if token != ('NA'):
                group = allele_to_group_map.get(token, None)
                new_tokens.append(token if group is None else group)
            else:
                new_tokens.append(token)
        print('\t'.join(new_tokens))
