#!/usr/bin/env python3
"""
Report Generator Module
Generates TXT and HTML reports from VCF files
"""

import os
from datetime import datetime


def parse_vcf(vcf_file):
    """Parse VCF file and return list of variants"""
    variants = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            variants.append({
                'chrom': parts[0],
                'pos': int(parts[1]),
                'ref': parts[3],
                'alt': parts[4],
                'qual': float(parts[5]),
                'type': 'INDEL' if 'INDEL' in parts[7] else 'SNP'
            })
    return variants


def generate_txt_report(variants, output_dir, sample_name):
    """Generate plain text report"""

    total = len(variants)
    snps = sum(1 for v in variants if v['type'] == 'SNP')
    indels = sum(1 for v in variants if v['type'] == 'INDEL')
    high_quality = sum(1 for v in variants if v['qual'] > 200)

    output_file = f"{output_dir}/{sample_name}_report.txt"

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write(f"  VARIANT DETECTION REPORT\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Sample Name:     {sample_name}\n")
        f.write(f"Reference:       Chromosome 17 (hg38)\n")
        f.write(f"Generated:       {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"VCF Source:      {sample_name}_variants.vcf\n\n")

        f.write("-" * 60 + "\n")
        f.write("  SUMMARY STATISTICS\n")
        f.write("-" * 60 + "\n")
        f.write(f"  Total Variants Detected:    {total}\n")
        f.write(f"    - SNPs:                    {snps}\n")
        f.write(f"    - INDELs:                  {indels}\n")
        f.write(f"    - High Quality (>200):     {high_quality}\n")
        f.write(f"    - Low Quality (≤200):      {total - high_quality}\n\n")

        # Quality distribution
        f.write("-" * 60 + "\n")
        f.write("  QUALITY DISTRIBUTION\n")
        f.write("-" * 60 + "\n")
        excellent = sum(1 for v in variants if v['qual'] >= 225)
        good = sum(1 for v in variants if 200 <= v['qual'] < 225)
        moderate = sum(1 for v in variants if 100 <= v['qual'] < 200)
        low = sum(1 for v in variants if v['qual'] < 100)

        f.write(f"  Excellent (≥225):   {excellent}\n")
        f.write(f"  Good (200-224):     {good}\n")
        f.write(f"  Moderate (100-199): {moderate}\n")
        f.write(f"  Low (<100):         {low}\n\n")

        # Top 15 variants
        f.write("-" * 60 + "\n")
        f.write("  TOP 15 VARIANTS\n")
        f.write("-" * 60 + "\n\n")

        sorted_variants = sorted(variants, key=lambda x: x['qual'], reverse=True)

        for i, var in enumerate(sorted_variants[:15]):
            f.write(f"  #{i + 1}\n")
            f.write(f"    Position:    {var['chrom']}:{var['pos']}\n")
            f.write(f"    Change:      {var['ref']} → {var['alt']}\n")
            f.write(f"    Quality:     {var['qual']:.2f}\n")
            f.write(f"    Type:        {var['type']}\n")
            f.write("\n")

        f.write("=" * 60 + "\n")
        f.write("  End of Report\n")
        f.write("=" * 60 + "\n")

    print(f"   TXT Report: {output_file}")
    return output_file


def generate_html_report(variants, output_dir, sample_name):
    """Generate HTML report"""

    total = len(variants)
    snps = sum(1 for v in variants if v['type'] == 'SNP')
    indels = sum(1 for v in variants if v['type'] == 'INDEL')
    high_quality = sum(1 for v in variants if v['qual'] > 200)

    output_file = f"{output_dir}/{sample_name}_report.html"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{sample_name} - Variant Report</title>
    <style>
        body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 30px; background: #f8f9fa; color: #2c3e50; }}
        .header {{ background: linear-gradient(135deg, #2c3e50, #3498db); color: white; padding: 25px; border-radius: 10px; margin-bottom: 25px; }}
        .header h1 {{ margin: 0; font-size: 26px; }}
        .header p {{ margin: 5px 0 0 0; opacity: 0.9; }}
        .stats {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 15px; margin-bottom: 25px; }}
        .stat-card {{ background: white; padding: 18px; border-radius: 8px; text-align: center; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        .stat-card .number {{ font-size: 30px; font-weight: bold; color: #3498db; }}
        .stat-card .label {{ font-size: 13px; color: #7f8c8d; margin-top: 5px; }}
        table {{ width: 100%; border-collapse: collapse; background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        th {{ background: #2c3e50; color: white; padding: 12px; text-align: left; }}
        td {{ padding: 10px 12px; border-bottom: 1px solid #ecf0f1; }}
        tr:hover {{ background: #f8f9fa; }}
        .footer {{ text-align: center; margin-top: 25px; color: #95a5a6; font-size: 13px; }}
        .badge {{ padding: 3px 8px; border-radius: 4px; font-size: 11px; font-weight: bold; }}
        .badge-snp {{ background: #3498db; color: white; }}
        .badge-indel {{ background: #e74c3c; color: white; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Variant Detection Report</h1>
        <p><strong>Sample:</strong> {sample_name} | <strong>Reference:</strong> Chromosome 17 (hg38) | <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>

    <div class="stats">
        <div class="stat-card">
            <div class="number">{total}</div>
            <div class="label">Total Variants</div>
        </div>
        <div class="stat-card">
            <div class="number">{snps}</div>
            <div class="label">SNPs</div>
        </div>
        <div class="stat-card">
            <div class="number">{indels}</div>
            <div class="label">INDELs</div>
        </div>
        <div class="stat-card">
            <div class="number">{high_quality}</div>
            <div class="label">High Quality (>200)</div>
        </div>
    </div>

    <h2>Top 20 Variants</h2>
    <table>
        <tr>
            <th>#</th>
            <th>Position</th>
            <th>Reference</th>
            <th>Variant</th>
            <th>Quality</th>
            <th>Type</th>
        </tr>
"""

    sorted_variants = sorted(variants, key=lambda x: x['qual'], reverse=True)

    for i, var in enumerate(sorted_variants[:20]):
        badge = 'badge-snp' if var['type'] == 'SNP' else 'badge-indel'
        html += f"""
        <tr>
            <td>{i + 1}</td>
            <td>{var['chrom']}:{var['pos']}</td>
            <td>{var['ref']}</td>
            <td><strong>{var['alt']}</strong></td>
            <td>{var['qual']:.2f}</td>
            <td><span class="badge {badge}">{var['type']}</span></td>
        </tr>"""

    html += """
    </table>

    <div class="footer">
        <p>Generated by DNA Sequencing Variant Detection Pipeline</p>
    </div>
</body>
</html>"""

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f"   HTML Report: {output_file}")
    return output_file


def generate_reports(vcf_file, output_dir, sample_name):
    """
    Generate all reports (TXT + HTML) from VCF file.

    Parameters:
    - vcf_file: Path to VCF file
    - output_dir: Directory to save reports
    - sample_name: Sample name for report titles
    """
    print(f"\n{'=' * 55}")
    print("GENERATING REPORTS")
    print(f"{'=' * 55}")

    variants = parse_vcf(vcf_file)

    print(f"\n{len(variants)} variants loaded from VCF\n")

    # Generate reports
    txt_file = generate_txt_report(variants, output_dir, sample_name)
    html_file = generate_html_report(variants, output_dir, sample_name)

    print(f"\nReports generated successfully!")

    return txt_file, html_file