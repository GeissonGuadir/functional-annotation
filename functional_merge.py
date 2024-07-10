import argparse
import csv
import os
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Process annotation data.')
    parser.add_argument('--filelist', required=True, help='File containing list of input files')
    parser.add_argument('--output_dir', required=True, help='Directory to save output TSV files')
    args = parser.parse_args()

    # Initialize dictionaries
    ncbi_dict = defaultdict(list)
    uniref_dict = defaultdict(list)
    signalp_dict = defaultdict(list)
    location_dict = defaultdict(dict)
    ipr_dict = defaultdict(set)
    go_dict = defaultdict(set)
    hmmer_dict = defaultdict(list)
    targetp_dict = defaultdict(str)
    ncbi_description_dict = {}

    # List for multiloc
    location_order = ['extracellular', 'ER', 'vacuolar', 'plasma membrane', 'Golgi apparatus', 'cytoplasmic', 'peroxisomal', 'mitochondrial', 'nuclear', 'chloroplast']

    # Location mapping
    location_mapping = {
        'extracellular': 'e',
        'ER': 'ER',
        'vacuolar': 'v',
        'plasma membrane': 'pm',
        'Golgi apparatus': 'Ga',
        'cytoplasmic': 'cy',
        'peroxisomal': 'p',
        'mitochondrial': 'm',
        'nuclear': 'n',
        'chloroplast': 'ch'
    }

    # Collect all unique genes
    all_genes = set()

    # Read the list file and categorize input files
    file_paths = {
        'ncbi': '',
        'uniref': '',
        'iprscan': '',
        'signalp': '',
        'targetp': '',
        'hmmer': '',
        'multiloc': '',
        'ncbidata': ''
    }

    with open(args.filelist, 'r') as filelist:
        for line in filelist:
            line = line.strip()
            if line.endswith('_ncbi_best_hit.out'):
                file_paths['ncbi'] = line
            elif line.endswith('_uniref_best_hit.out'):
                file_paths['uniref'] = line
            elif line.endswith('.iprscan.tsv'):
                file_paths['iprscan'] = line
            elif line.endswith('.gff3'):
                file_paths['signalp'] = line
            elif line.endswith('_summary.targetp2'):
                file_paths['targetp'] = line
            elif line.endswith('.domtblout'):
                file_paths['hmmer'] = line
            elif line.endswith('_ml2.txt'):
                file_paths['multiloc'] = line
            elif line.endswith('.prot'):
                file_paths['ncbidata'] = line

    # NCBI description
    print("Processing NCBI descriptions...")
    with open(file_paths['ncbidata'], 'r') as ncbi_file:
        for line in ncbi_file:
            if line.startswith('>'):  # Check if the line starts with '>'
                parts = line[1:].strip().split(' ', 1)  # Skip the '>' and split into two parts
                accession = parts[0]  # The accession number is the first part
                description = parts[1] if len(parts) > 1 else "Description not found"
                ncbi_description_dict[accession] = description

    # Multiloc
    print("Processing MultiLoc data...")
    with open(file_paths['multiloc'], 'r') as multiloc:
        for line in multiloc:
            gene, *locations = line.strip().split('\t')
            all_genes.add(gene)
            for location in locations:
                if ':' in location:
                    name, value = location.split(':')
                    location_dict[gene][name.strip()] = value.strip()

    # SignalP
    print("Processing SignalP data...")
    with open(file_paths['signalp'], 'r') as signalp:
        for line in signalp:
            if not line.startswith('#'):
                fields = line.strip().split()
                gene = fields[0]
                signalp_dict[gene] = fields[4:6]
                all_genes.add(gene)

    # TargetP
    print("Processing TargetP data...")
    with open(file_paths['targetp'], 'r') as targetpfile:
        for line in targetpfile:
            if not line.startswith('#'):
                fields = line.strip().split()
                gene = fields[0]
                targetp_dict[gene] = " ".join(fields[1:])
                all_genes.add(gene)

    # ncbi and uniref
    print("Processing NCBI and Uniref data...")
    for file_path in [file_paths['ncbi'], file_paths['uniref']]:
        with open(file_path, 'r') as file:
            for line in file:
                if not line.startswith('#'):
                    fields = line.strip().split()
                    gene = fields[0]
                    all_genes.add(gene)
                    if file_path == file_paths['ncbi']:
                        ncbi_dict[gene] = [fields[1], fields[11]]
                    else:
                        uniref_dict[gene] = [fields[1], fields[11]]

    # hmmer
    print("Processing HMMER data...")
    with open(file_paths['hmmer'], 'r') as hmmerfile:
        for line in hmmerfile:
            if not line.startswith('#'):
                fields = line.strip().split()
                gene = fields[3]
                all_genes.add(gene)
                accession = fields[1]
                hmmer_dict[gene].append(accession)

    # iprscan
    print("Processing IPRScan data...")
    with open(file_paths['iprscan'], 'r') as iprscanfile:
        for line in iprscanfile:
            fields = line.strip().split('\t')
            gene = fields[0]
            all_genes.add(gene)
            ipr = fields[11] if fields[11] != '-' else ''
            go_terms = fields[13]
            if ipr:
                ipr_dict[gene].add(ipr)
            if go_terms:
                for go_term in go_terms.split('|'):
                    if 'GO:' in go_term:
                        go = go_term
                        go_dict[gene].add(go)


    print("Total genes collected:", len(all_genes))

    # Determine output file name
    output_file = os.path.join(args.output_dir, os.path.basename(file_paths['ncbi']).replace('_ncbi_best_hit.out', '.tsv'))
    print("Writing results to:", output_file)
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')

        header = ["ID sequence", "NCBI DESCRIPTION", "NCBI Subject; Bit score", "Uniref Subject; Bit score", "SignalP Pos; Pr", "TargetP Prediction; noTP; SP; mTP; cTP; luTP; CS Position",
                  "; ".join([location_mapping[name] for name in location_order]), "IPRSCAN GO", "IPRSCAN IPR", "Hmmer Pfam"]
        writer.writerow(header)

        for gene in all_genes:
            results = [location_dict[gene].get(name, "N/A") for name in location_order]

            ncbi_values = ncbi_dict.get(gene, ["N/A"])
            description = ncbi_description_dict.get(ncbi_values[0].split(';')[0], "Description not found")

            targetp_info = targetp_dict[gene].replace(' ', ';') if gene in targetp_dict else "N/A"
            ipr_terms = ";".join(sorted(ipr_dict[gene]))
            go_terms = ";".join(sorted(go_dict[gene]))

            row = [gene, description, ";".join(ncbi_values), ";".join(uniref_dict.get(gene, ["N/A"])),
                   ";".join(signalp_dict.get(gene, ["N/A"])), targetp_info, ";".join(results), go_terms, ipr_terms, ";".join(hmmer_dict.get(gene, ["N/A"]))]
            writer.writerow(row)

if __name__ == '__main__':
    main()
