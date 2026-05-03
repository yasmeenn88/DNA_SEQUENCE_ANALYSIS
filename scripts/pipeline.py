#!/usr/bin/env python3
"""
Variant Detection Pipeline - Core Module
Runs: minimap2 -> SAMtools -> BCFtools
"""

import subprocess
import sys


def run_command(cmd, step_name):
    """Run a shell command and exit if it fails"""
    print(f"\n{'=' * 55}")
    print(f"[{step_name}]")
    print(f"{'=' * 55}")

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.stdout:
        print(result.stdout.strip())

    if result.returncode != 0:
        print(f"ERROR in: {step_name}")
        if result.stderr:
            print(result.stderr.strip())
        sys.exit(1)

    print(f"DONE: {step_name}")


def run_pipeline(reference, fastq1, fastq2, output_prefix, vcf_output):
    """
    Run complete variant detection pipeline.

    Parameters:
    - reference: Path to reference genome (FASTA)
    - fastq1: Path to first FASTQ file
    - fastq2: Path to second FASTQ file
    - output_prefix: Prefix for intermediate files
    - vcf_output: Final VCF output path
    """

    print("\n" + "=" * 55)
    print("VARIANT DETECTION PIPELINE")
    print("=" * 55)
    print(f"Reference: {reference}")
    print(f"Reads: {fastq1} + {fastq2}")
    print(f"Output: {vcf_output}")
    print("=" * 55)

    # Step 1: Index reference
    run_command(f"minimap2 -d {reference}.mmi {reference}", "1/6 - Indexing Reference Genome")

    # Step 2: Alignment
    run_command(
        f"minimap2 -ax sr {reference} {fastq1} {fastq2} > {output_prefix}_aligned.sam",
        "2/6 - Sequence Alignment (minimap2)"
    )

    # Step 3: SAM to BAM
    run_command(
        f"samtools view -bS {output_prefix}_aligned.sam > {output_prefix}_aligned.bam",
        "3/6 - Converting SAM to BAM"
    )

    # Step 4: Sort BAM
    run_command(
        f"samtools sort {output_prefix}_aligned.bam -o {output_prefix}_sorted.bam",
        "4/6 - Sorting BAM"
    )

    # Step 5: Index BAM
    run_command(
        f"samtools index {output_prefix}_sorted.bam",
        "5/6 - Indexing BAM"
    )

    # Step 6: Variant Calling
    run_command(
        f"bcftools mpileup -f {reference} {output_prefix}_sorted.bam | bcftools call -mv -o {vcf_output}",
        "6/6 - Variant Calling (BCFtools)"
    )

    # Quick summary
    count_cmd = f"grep -v '^#' {vcf_output} | wc -l"

    result = subprocess.run(count_cmd, shell=True, capture_output=True, text=True)
    variant_count = result.stdout.strip()

    print(f"\n{'=' * 55}")
    print(f"PIPELINE COMPLETE!")
    print(f"Total Variants Found: {variant_count}")
    print(f"VCF File: {vcf_output}")
    print(f"{'=' * 55}\n")
