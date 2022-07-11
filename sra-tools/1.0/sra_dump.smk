import pandas as pd
import os

configfile: "test_sra_dump.yaml"

accessions = pd.read_csv(config["accessions"], names = ["accession"])

in_file_type = config["in_file_type"]
assert in_file_type in ["sam", "fastq"], (
    f"Unrecognized in_file_type {in_file_type}. \n"
    "Please specify one of [\"sam\", \"fastq\"]. "
)

wildcard_constraints:
    accession = "|".join(accessions["accession"].tolist())

def get_ngc_arg(ngc = config["params"]["ngc"]):
    if ngc not in ["__UPDATE__", ""]:
        assert os.path.exists(ngc), (
            f"ngc file {ngc} does not exist. \n"
            "Please specify a .ngc file or leave config[\"params\"][\"ngc\"] blank if the repository is public. "
        )
        ngc_param = f"--ngc {ngc}"
    else:
        ngc_param = ""
    return ngc_param

checkpoint fasterq_dump:
    output:
        complete = touch("{accession}.fastq.complete"),
    params:
        ngc = get_ngc_arg(),
        opts = config["params"]["fastq_dump"]
    conda:
        config["conda"]["sra-tools"]
    threads: config["threads"]["fastq_dump"]
    resources:
        **config["resources"]["fastq_dump"]
    wildcard_constraints:
        accession = "|".join(accessions["accession"].tolist())
    group: "fasterq"
    shell:
        "fasterq-dump -e {threads} {params.opts} {params.ngc} {wildcards.accession}"

rule gzip_fastq:
    input:
        fastq = "{accession}{read}.fastq"
    output:
        fastq = "{accession}{read}.fastq.gz"
    threads: 1
    resources:
        **config["resources"]["fastq_dump"]
    wildcard_constraints:
        accession = "|".join(accessions["accession"].tolist())
    group: "fasterq"
    shell:
        "gzip -c {input.fastq} > {output.fastq} && rm -f {input.fastq}"

def get_fastq_read(wildcards):
    # Get the path to the output directory from the checkpoint fasterq_dump
    # **wildcards unpacks the wildcards object from the rule where this input function is called
    checkpoint_output = checkpoints.fasterq_dump.get(accession=wildcards.accession).output.complete
    print(glob_wildcards(checkpoint_output.replace(".fastq.complete", "") + "{read}.fastq"))
    # Obtain the read wildcards using the glob_wildcards function
    fastqs = expand(
        rules.gzip_fastq.output.fastq,
        accession=wildcards.accession,
        read = glob_wildcards(checkpoint_output.replace(".fastq.complete", "") + "{read}.fastq").read
        )
    return fastqs

rule fastq_complete:
    input:
        get_fastq_read
    output:
        touch("{accession}.gzip.complete")

def get_genome_arg(
    out_type = config["out_file_type"],
    ref_fasta = config["params"]["ref"]
):
    if out_type == "cram":
        assert os.path.exists(ref_fasta), (f"Reference fasta file {ref_fasta} does not exist. ")
        param = f"-CS -T {ref_fasta}"
    elif out_type == "bam":
        param = "-bS"
    else:
        raise AssertionError(f"Specified out_file_type is \"{out_type}\", not one of [\"bam\", \"cram\"]. ")
    return(param)


rule sam_dump:
    output:
        sam = pipe("{accession}.sam")
    params:
        ngc = get_ngc_arg(),
        opts = config["params"]["sam-dump"]
    conda:
        config["conda"]["sra-tools"]
    threads: 1
    resources:
        mem_mb = 10000
    shell:
        "sam-dump {params.opts} {wildcards.accession} > {output.sam}"

rule sam_to_bam:
    input:
        sam = str(rules.sam_dump.output.sam)
    output:
        bam = "{accession}.hg19." + config["out_file_type"],
        complete = touch("{accession}.complete")
    params:
        bam_cram = get_genome_arg(),
        opts = config["params"]["samtools"]
    conda:
        config["conda"]["sra-tools"]
    threads: config["threads"]["sam_to_bam"]
    resources:
        **config["resources"]["sam_to_bam"]
    shell:
        "samtools view {params.bam_cram} -@ {threads} {input.sam} > {output.bam}"

rule index_bam:
    input:
        bam = str(rules.sam_to_bam.output.bam)
    output:
        complete = touch("{accession}.indexed")
    conda:
        config["conda"]["sra-tools"]
    threads: config["threads"]["index_bam"]
    resources:
        **config["resources"]["index_bam"]
    shell:
        "samtools index -@ {threads} {input.bam}"


rule all:
    input:
        expand(
            [
                str(rules.sam_to_bam.output.bam),
                str(rules.sam_to_bam.output.complete),
                str(rules.index_bam.output.complete)
            ],
            accession = accessions["accession"]
        ) if in_file_type == "sam" else expand(
            str(rules.fastq_complete.output[0]),
            accession = accessions["accession"]
        )