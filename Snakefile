configfile: "config.yaml"

import os

if '.' in config['fasta_file']:
    species_name = config['fasta_file'][:config['fasta_file'].rindex('.')]
else:
    species_name = config['fasta_file']

# Define the final output expected at the end of the pipeline
rule all:
    input:
        expand(f"functional_outputs/{species_name}.tsv")

rule agat:
    input:
        final_filtering_gff=f"inputs/{config['final_filtering_gff']}"
    output:
        agat_output=f"agat_outputs/{species_name}.gff"
    params:
        volume_name = config["volume_name"]
    singularity:
        "images/agat.sif"
    shell:
        "agat_sp_keep_longest_isoform.pl -gff {params.volume_name}/{input} -o {output}"

rule gff_read:
    input:
        agat_output=rules.agat.output,
        fasta_file=f"inputs/{config['fasta_file']}"
    output:
        aa=f"gff_output/{species_name}.aa",
        cds=f"gff_output/{species_name}.cds"
    params:
        volume_name=config['volume_name']
    singularity:
        "images/gffread.sif"
    shell:
        """
        gffread {params.volume_name}/{input.agat_output} -g {params.volume_name}/{input.fasta_file} -J -S -y {params.volume_name}/{output.aa}
        gffread {params.volume_name}/{input.agat_output} -g {params.volume_name}/{input.fasta_file} -J -S -x {params.volume_name}/{output.cds}
        """

rule mmseqs2:
    input:
        gff_output=rules.gff_read.output,
        aa=f"gff_output/{species_name}.aa"
    output:
         ncbi=f"mmseqs_output/{species_name}_ncbi_select.ofmt6",
         uniref=f"mmseqs_output/{species_name}_uniref50.ofmt6",
    params:
         volume_name=config['volume_name']
    singularity:
        "images/mmseqs2.sif"
    shell:
        """
        mmseqs easy-search {params.volume_name}/{input.aa} {params.volume_name}/databases/current/GenomesDatabase.v1.prot.mmseqs {output.ncbi} /tmp --threads 20 --format-mode 0 --start-sens 2 -s 7 --sens-steps 3
        mmseqs easy-search {params.volume_name}/{input.aa} {params.volume_name}/databases/mmseqs_db/uniref50 {output.uniref} /tmp --threads 20 --format-mode 0 --start-sens 2 -s 7 --sens-steps 3
        """

rule best_hits:
    input:
        mmseqs2=rules.mmseqs2.output,
        ncbi=f"mmseqs_output/{species_name}_ncbi_select.ofmt6",
        uniref=f"mmseqs_output/{species_name}_uniref50.ofmt6"
    output:
        ncbi_best_hit=f"mmseqs_output/{species_name}_ncbi_best_hit.out",
        uniref_best_hit=f"mmseqs_output/{species_name}_uniref_best_hit.out"
    params:
        volume_name=config['volume_name']
    shell:
        """
        sort -k12 -t $'\t' -nr {input.ncbi} | awk -F "\t" ' ! a[$1]++ && $11 < 1e-5' > {output.ncbi_best_hit}
        sort -k12 -t $'\t' -nr {input.uniref} | awk -F "\t" ' ! a[$1]++ && $11 < 1e-5' > {output.uniref_best_hit}
        """

rule interproscan:
    input:
        gff_output=rules.best_hits.output,
        aa=f"gff_output/{species_name}.aa"
    output:
        iprscan_output=f"iprscan_output/{species_name}.iprscan.tsv"
    params:
        volume_name=config['volume_name']
    shell:
        """
        singularity exec -B {params.volume_name}/databases/interproscan-5.67-99.0/data:/opt/interproscan/data -B {params.volume_name}:/data {params.volume_name}/images/iprscan.sif /opt/interproscan/interproscan.sh --formats GFF3 TSV  --goterms --pathways --iprlookup --input /data/{input.aa} --cpu 30 --output-file-base /data/iprscan_output/{species_name}.iprscan
        """

rule prot_list:
    input:
        inteproscan_output=rules.interproscan.output,
        aa=f"gff_output/{species_name}.aa",
        iprscan_output=f"iprscan_output/{species_name}.iprscan.tsv"
    output:
        multigo=f"multiloc_outputs/{species_name}_multilocGo.out",
        prot_list=f"multiloc_outputs/{species_name}_prot.list"
    params:
        volume_name=config['volume_name']
    shell:
        """
        awk -F"\t" '$14 ~ /GO/ {{print $1,$14}}' {params.volume_name}/{input.iprscan_output} | sed 's/|/ /g' > {output.multigo}
        egrep \> {params.volume_name}/{input.aa} | cut -b 2- | awk '{{print $1}}' > {output.prot_list}
        mkdir splits
        cd splits
        split -n l/20 {params.volume_name}/{output.prot_list} #20 is the number of threads
        """

rule multiloc:
    input:
        prot_list_output=rules.prot_list.output,
        multiloc_script=f"scripts/{config['multiloc_script']}",
        aa=f"gff_output/{species_name}.aa",
        multigo=f"multiloc_outputs/{species_name}_multilocGo.out"
    output:
        multiloc_output=f"multiloc_outputs/{species_name}_ml2.txt"
    params:
        volume_name=config['volume_name']
    shell:
        """
        cd splits
        declare -a pid_array
        for x in x{{a..z}}{{a..z}}; do
            mkdir $x"_dir";
            {params.volume_name}/{input.multiloc_script} $x {input.aa} {input.multigo} $x"_dir" {params.volume_name} & pid_array+=("$!");
        done
            for pid in "${{pid_array[@]}}"; do
            wait $pid
        done
        cat x*/ml2.out > {params.volume_name}/{output.multiloc_output}
        cd -
        rm -rf splits        
        """

rule hmmer:
    input:
        gff_output=rules.multiloc.output,
        aa=f"gff_output/{species_name}.aa"
    output:
        hmmer_output=f"hmmer_outputs/{species_name}.domtblout"
    params:
        volume_name=config['volume_name']
    singularity:
        "images/hmmer3.sif"
    shell:
        """
        hmmscan --cpu 4 --domtblout {params.volume_name}/{output.hmmer_output} {params.volume_name}/databases/Pfam-A.hmm {params.volume_name}/{input.aa}
        """

rule sigtarp:
    input:
        gff_output=rules.hmmer.output,
        aa=f"gff_output/{species_name}.aa"
    output:
        signalp_output=f"sigtarp_outputs/{species_name}_summary.signalp5",
        targetp_output=f"sigtarp_outputs/{species_name}_summary.targetp2"
    params:
        volume_name=config['volume_name']
    singularity:
        "images/sigtarp.sif"
    shell:
        """
        cd sigtarp_outputs
        signalp -fasta {params.volume_name}/{input.aa} -org euk -format short -gff3
        targetp -fasta {params.volume_name}/{input.aa} -org pl -format short
        cd -
        """

rule generate_tsv:
    input:
	    sigtar=rules.sigtarp.output,
        mmseqs=rules.best_hits.output,
        iprs=rules.interproscan.output,
        hmmer=rules.hmmer.output,
        multiloc=rules.multiloc.output,
        ncbi_database=f"databases/current/GenomesDatabase.v1.prot"
    output:
	    output_file=f"functional_outputs/Erigeron_canadensis.tsv",
        outputs_list=f"functional_outputs/outputs_list"
    params:
	    volume_name=config['volume_name']
    shell:
        """
        for item in {input.mmseqs} {input.iprs} {input.hmmer} {input.multiloc} {input.sigtar} {input.ncbi_database}
        do
            echo {params.volume_name}/$item >> {params.volume_name}/{output.outputs_list}
        done
        python {params.volume_name}/scripts/functional_merge.py --filelist {params.volume_name}/{output.outputs_list} --output_dir {params.volume_name}/functional_outputs
        """
