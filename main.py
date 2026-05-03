#!/usr/bin/env python3
"""
🧬 DNA Sequencing Variant Detection Pipeline
Main Entry Point - Run with: python3 main.py --sample <sample_name>
"""

import sys
import os
import argparse

# Add scripts folder to path
sys.path.append('scripts')

from scripts.pipeline import run_pipeline
from scripts.clinvar_api import validate_with_clinvar
from scripts.report_generator import generate_reports
from scripts.plot_generator import generate_plots


def create_sample_dirs(sample_name):
    """Create all directories for a sample"""
    dirs = [
        f"samples/{sample_name}",
        f"samples/{sample_name}/raw",
        f"samples/{sample_name}/aligned",
        f"samples/{sample_name}/variants",
        f"samples/{sample_name}/reports",
        f"samples/{sample_name}/plots",
        f"samples/{sample_name}/validation"
    ]
    for d in dirs:
        os.makedirs(d, exist_ok=True)
    print(f"📁 Created directories for: {sample_name}")


def find_fastq_files(sample_name):
    """Auto-detect FASTQ files in sample/raw/ folder"""
    raw_dir = f"samples/{sample_name}/raw"

    if not os.path.exists(raw_dir):
        print(f"❌ Directory not found: {raw_dir}")
        sys.exit(1)

    files = sorted(os.listdir(raw_dir))
    fastq_files = [f for f in files if f.endswith('.fastq') or f.endswith('.fq')]

    if len(fastq_files) >= 2:
        fastq1 = f"{raw_dir}/{fastq_files[0]}"
        fastq2 = f"{raw_dir}/{fastq_files[1]}"
        print(f"📁 Auto-detected FASTQ files:")
        print(f"   R1: {fastq1}")
        print(f"   R2: {fastq2}")
        return fastq1, fastq2
    elif len(fastq_files) == 1:
        fastq1 = f"{raw_dir}/{fastq_files[0]}"
        print(f"⚠️  Only 1 FASTQ file found (single-end mode)")
        return fastq1, None
    else:
        print(f"❌ No FASTQ files found in {raw_dir}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='🧬 DNA Sequencing Variant Detection Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 main.py --sample brca1
  python3 main.py --sample tp53 --reference reference/chr17.fa
  python3 main.py --sample brca1 --fastq1 data_R1.fastq --fastq2 data_R2.fastq
        """
    )

    parser.add_argument('--sample', '-s', required=True,
                        help='Sample name (e.g., brca1, tp53)')
    parser.add_argument('--reference', '-r', default='reference/chr17.fa',
                        help='Path to reference genome (default: reference/chr17.fa)')
    parser.add_argument('--fastq1', '-f1', default=None,
                        help='Path to first FASTQ file (auto-detected if not provided)')
    parser.add_argument('--fastq2', '-f2', default=None,
                        help='Path to second FASTQ file (auto-detected if not provided)')
    parser.add_argument('--skip-clinvar', action='store_true',
                        help='Skip ClinVar validation')
    parser.add_argument('--skip-plots', action='store_true',
                        help='Skip plot generation')

    args = parser.parse_args()

    sample_name = args.sample
    reference = args.reference

    # ==========================================
    # STEP 1: Setup directories
    # ==========================================
    print("\n" + "=" * 55)
    print("🧬 DNA SEQUENCING VARIANT DETECTION PIPELINE")
    print("=" * 55)
    print(f"Sample:    {sample_name}")
    print(f"Reference: {reference}")
    print("=" * 55)

    create_sample_dirs(sample_name)

    # ==========================================
    # STEP 2: Find FASTQ files
    # ==========================================
    if args.fastq1:
        fastq1 = args.fastq1
        fastq2 = args.fastq2 or fastq1  # If no fastq2, use fastq1
    else:
        fastq1, fastq2 = find_fastq_files(sample_name)

    if not os.path.exists(fastq1):
        print(f"❌ FASTQ file not found: {fastq1}")
        sys.exit(1)

    if fastq2 and not os.path.exists(fastq2):
        print(f"❌ FASTQ file not found: {fastq2}")
        sys.exit(1)

    # ==========================================
    # STEP 3: Define output paths
    # ==========================================
    output_prefix = f"samples/{sample_name}/aligned/{sample_name}"
    vcf_output = f"samples/{sample_name}/variants/{sample_name}_variants.vcf"
    reports_dir = f"samples/{sample_name}/reports"
    plots_dir = f"samples/{sample_name}/plots"
    validation_dir = f"samples/{sample_name}/validation"

    # ==========================================
    # STEP 4: Run Pipeline
    # ==========================================
    if fastq2:
        run_pipeline(reference, fastq1, fastq2, output_prefix, vcf_output)
    else:
        print("⚠️ Single-end mode: treating as single FASTQ")
        # For single-end, duplicate fastq1 as both inputs
        run_pipeline(reference, fastq1, fastq1, output_prefix, vcf_output)

    # ==========================================
    # STEP 5: ClinVar Validation
    # ==========================================
    if not args.skip_clinvar:
        validate_with_clinvar(vcf_output, validation_dir)
    else:
        print("\n⏭️  Skipping ClinVar validation")

    # ==========================================
    # STEP 6: Generate Reports
    # ==========================================
    generate_reports(vcf_output, reports_dir, sample_name)

    # ==========================================
    # STEP 7: Generate Plots
    # ==========================================
    if not args.skip_plots:
        generate_plots(vcf_output, plots_dir, sample_name)
    else:
        print("\n⏭️  Skipping plot generation")

    # ==========================================
    # FINAL SUMMARY
    # ==========================================
    print("\n" + "=" * 55)
    print("🎉 PIPELINE COMPLETE!")
    print("=" * 55)
    print(f"""
📁 Results for: {sample_name}
   ├── aligned/       → SAM + BAM files
   ├── variants/      → {sample_name}_variants.vcf
   ├── reports/       → TXT + HTML reports
   ├── plots/         → Manhattan + Heatmap + Circos
   └── validation/    → clinvar_validation.json

✅ All done! Check: samples/{sample_name}/
""")
    print("=" * 55)


if __name__ == "__main__":
    main()



















# # import pandas as pd
# #
# #
# # def parse_vcf(vcf_file):
# #     """Parse VCF file and return list of variants"""
# #     variants = []
# #     with open(vcf_file, 'r') as f:
# #         for line in f:
# #             if line.startswith('#'):
# #                 continue
# #             parts = line.strip().split('\t')
# #             variant = {
# #                 'CHROM': parts[0],
# #                 'POS': int(parts[1]),
# #                 'REF': parts[3],
# #                 'ALT': parts[4],
# #                 'QUAL': float(parts[5]),
# #                 'TYPE': 'INDEL' if 'INDEL' in parts[7] else 'SNP'
# #             }
# #             variants.append(variant)
# #     return variants
# #
# #
# # def generate_report(vcf_file, output_file):
# #     """Generate summary report from VCF file"""
# #     variants = parse_vcf(vcf_file)
# #     df = pd.DataFrame(variants)
# #
# #     # Statistics
# #     total = len(df)
# #     snps = len(df[df['TYPE'] == 'SNP'])
# #     indels = len(df[df['TYPE'] == 'INDEL'])
# #     high_quality = len(df[df['QUAL'] > 200])
# #
# #     # Write report
# #     with open(output_file, 'w') as f:
# #         f.write("=" * 60 + "\n")
# #         f.write("BRCA1 MUTATION DETECTION REPORT\n")
# #         f.write("=" * 60 + "\n\n")
# #         f.write(f"Reference Genome: Chromosome 17 (hg38)\n\n")
# #         f.write("SUMMARY STATISTICS:\n")
# #         f.write("-" * 40 + "\n")
# #         f.write(f"Total Variants Detected: {total}\n")
# #         f.write(f"  - Single Nucleotide Polymorphisms (SNPs): {snps}\n")
# #         f.write(f"  - Insertions/Deletions (INDELs): {indels}\n")
# #         f.write(f"  - High Quality (QUAL > 200): {high_quality}\n\n")
# #
# #         # Top 10 variants
# #         f.write("TOP 10 VARIANTS:\n")
# #         f.write("-" * 40 + "\n")
# #         for _, var in df.head(10).iterrows():
# #             f.write(f"Position: {var['POS']}\n")
# #             f.write(f"  Reference: {var['REF']}\n")
# #             f.write(f"  Variant: {var['ALT']}\n")
# #             f.write(f"  Quality: {var['QUAL']:.2f}\n")
# #             f.write(f"  Type: {var['TYPE']}\n\n")
# #
# #     print(f"Report generated: {output_file}")
# #     print(f"Total variants: {total}")
# #
# #
# # if __name__ == "__main__":
# #     generate_report("brca1_variants.vcf", "mutation_report.txt")
#
#
# python3 << 'EOF'
# import requests
# import time
#
# variants = [
#     ("chr17", 7563481, "C", "G"),
#     ("chr17", 7563997, "G", "C"),
#     ("chr17", 7564955, "A", "C"),
#     ("chr17", 7568207, "C", "G"),
# ]
#
# print("🔍 Searching ClinVar for TP53 mutations...\n")
#
# for chrom, pos, ref, alt in variants:
#     search_term = f"{chrom}:{pos}{ref}>{alt}"
#     url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
#     params = {
#         'db': 'clinvar',
#         'term': search_term,
#         'retmax': 1,
#         'retmode': 'json'
#     }
#
#     try:
#         response = requests.get(url, params=params, timeout=10)
#         data = response.json()
#         id_list = data.get('esearchresult', {}).get('idlist', [])
#
#         if id_list:
#             print(f"✅ {search_term}: FOUND in ClinVar (ID: {id_list[0]})")
#         else:
#             print(f"❌ {search_term}: Not in ClinVar")
#     except Exception as e:
#         print(f"⚠️ {search_term}: Error - {e}")
#
#     time.sleep(0.5)
#
# print("\n✅ Done!")
# EOF