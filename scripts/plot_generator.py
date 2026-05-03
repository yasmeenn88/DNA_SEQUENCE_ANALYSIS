"""
Plot Generator Module
Generates: Manhattan Plot + Heatmap + Circos Plot
"""

import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import matplotlib

matplotlib.use('Agg')  # Use non-interactive backend

plt.rcParams['font.family'] = 'DejaVu Sans'


def parse_vcf(vcf_file):
    """Parse VCF file and return lists of positions, qualities, types"""
    positions = []
    qualities = []
    types = []
    ref_bases = []
    alt_bases = []

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            qual = float(parts[5])
            var_type = 'INDEL' if 'INDEL' in parts[7] else 'SNP'

            positions.append(pos)
            qualities.append(qual)
            types.append(var_type)
            ref_bases.append(ref)
            alt_bases.append(alt)

    return positions, qualities, types, ref_bases, alt_bases


def manhattan_plot(vcf_file, output_dir, sample_name):
    """Generate Manhattan plot
    - x-axis -> position on chr17.
    - y-axis -> quality score
    النقط العالية بتكون بنسبة كبيرة حقيقية مش نويز
    ولو كلاستر متجمع عنده نقط كتير بعتبره هوت سبوت

    """
    positions, qualities, types, _, _ = parse_vcf(vcf_file)

    if not positions:
        print("   No variants found for Manhattan plot")
        return

    fig, ax = plt.subplots(figsize=(14, 6))

    # Separate SNPs and INDELs
    snp_pos = [p for p, t in zip(positions, types) if t == 'SNP']
    snp_qual = [q for q, t in zip(qualities, types) if t == 'SNP']
    indel_pos = [p for p, t in zip(positions, types) if t == 'INDEL']
    indel_qual = [q for q, t in zip(qualities, types) if t == 'INDEL']

    # Plot
    ax.scatter(snp_pos, snp_qual, c='#3498db', alpha=0.7, s=60,
               edgecolors='white', linewidth=0.5, label=f'SNPs ({len(snp_pos)})')
    ax.scatter(indel_pos, indel_qual, c='#e74c3c', alpha=0.8, s=80,
               edgecolors='white', linewidth=0.5, marker='^', label=f'INDELs ({len(indel_pos)})')

    # Threshold lines
    ax.axhline(y=225, color='#f39c12', linestyle='--', alpha=0.7, label='High (225)')
    ax.axhline(y=200, color='#95a5a6', linestyle=':', alpha=0.5, label='Medium (200)')

    ax.set_xlabel('Position on Chromosome 17 (Mb)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Variant Quality Score (Phred)', fontsize=13, fontweight='bold')
    ax.set_title(f'Manhattan Plot - {sample_name}\nChromosome 17 (hg38)',
                 fontsize=16, fontweight='bold', pad=15)

    # Format x-axis to Mb
    ticks = ax.get_xticks()
    ax.set_xticklabels([f'{x / 1e6:.1f}' for x in ticks])

    ax.legend(loc='upper right', fontsize=11, framealpha=0.9)
    ax.grid(alpha=0.3, linestyle='--')

    plt.tight_layout()
    output_file = f"{output_dir}/manhattan_plot.png"
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close()

    print(f"   Manhattan Plot: {output_file}")


# ال heatmap متقسمة ل Histogram وDensity + خطوط و Base Substitution
def mutation_heatmap(vcf_file, output_dir, sample_name):
    """Generate heatmap-style plot"""
    positions, qualities, types, ref_bases, alt_bases = parse_vcf(vcf_file)

    if not positions:
        print("   No variants found for Heatmap")
        return

    fig, axes = plt.subplots(3, 1, figsize=(14, 10),
                             gridspec_kw={'height_ratios': [1, 2, 1]})

    # Plot 1: Density histogram
    ax1 = axes[0]
    ax1.hist(positions, bins=50, color='#3498db', alpha=0.7, edgecolor='white')
    ax1.set_ylabel('Density', fontsize=12, fontweight='bold')
    ax1.set_title(f'Mutation Distribution & Density - {sample_name}\nChromosome 17',
                  fontsize=16, fontweight='bold')
    ax1.grid(alpha=0.3)
    ax1.set_xticklabels([])

    # Plot 2: Strip plot with density
    # خيوط كتير-> منطقة غنية ميوتيشن
    ax2 = axes[1]
    hist, bin_edges = np.histogram(positions, bins=80)
    max_density = hist.max() if hist.max() > 0 else 1

    for pos, var_type in zip(positions, types):
        bin_idx = np.digitize(pos, bin_edges) - 1
        bin_idx = max(0, min(bin_idx, len(hist) - 1))
        density = hist[bin_idx]
        alpha_val = 0.2 + (density / max_density) * 0.8

        color = '#e74c3c' if var_type == 'SNP' else '#2ecc71'
        linewidth = 1.5 if var_type == 'SNP' else 2.5
        ax2.axvline(x=pos, color=color, alpha=alpha_val, linewidth=linewidth)

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    ax2_twin = ax2.twinx()
    ax2_twin.plot(bin_centers, hist, color='#3498db', linewidth=2, alpha=0.5)
    ax2_twin.fill_between(bin_centers, 0, hist, color='#3498db', alpha=0.1)
    ax2_twin.set_ylabel('Density', fontsize=11, color='#3498db')

    ax2.set_yticklabels([])
    ax2.grid(alpha=0.2, axis='x')

    from matplotlib.patches import Patch
    ax2.legend(handles=[
        Patch(facecolor='#e74c3c', alpha=0.7, label=f'SNPs ({types.count("SNP")})'),
        Patch(facecolor='#2ecc71', alpha=0.7, label=f'INDELs ({types.count("INDEL")})')
    ], loc='upper right')

    # Plot 3: Base change frequencies
    # لو نوع من التغيرات السينجل متكرر جدا بتبقا دلالة ع مرض معين
    ax3 = axes[2]
    base_changes = Counter()
    for i, t in enumerate(types):
        if t == 'SNP' and len(ref_bases[i]) == 1 and len(alt_bases[i]) == 1:
            base_changes[f"{ref_bases[i]}→{alt_bases[i]}"] += 1

    if base_changes:
        changes = list(base_changes.keys())[:10]
        counts = [base_changes[c] for c in changes]
        colors_palette = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6',
                          '#1abc9c', '#e67e22', '#34495e', '#95a5a6', '#c0392b']

        bars = ax3.bar(range(len(changes)), counts, color=colors_palette[:len(changes)],
                       edgecolor='white', linewidth=2)
        ax3.set_xticks(range(len(changes)))
        ax3.set_xticklabels(changes, fontsize=11, fontweight='bold')
        ax3.set_ylabel('Count', fontsize=12, fontweight='bold')
        ax3.set_xlabel('Base Substitution', fontsize=12, fontweight='bold')
        ax3.set_title('Most Common Base Substitutions', fontsize=13, fontweight='bold')
        ax3.grid(alpha=0.2, axis='y')

        for bar, count in zip(bars, counts):
            ax3.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15,
                     str(count), ha='center', fontweight='bold', fontsize=10)

    plt.tight_layout()
    output_file = f"{output_dir}/heatmap.png"
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close()

    print(f"   Heatmap: {output_file}")


def circos_plot(vcf_file, output_dir, sample_name):
    """Generate circular plot"""
    positions, qualities, types, _, _ = parse_vcf(vcf_file)

    if not positions:
        print("   No variants found for Circos plot")
        return

    chr17_size = 83257441  # Chromosome 17 size in bp

    fig, ax = plt.subplots(figsize=(12, 12), subplot_kw={'projection': 'polar'})

    for pos, qual, var_type in zip(positions, qualities, types):
        angle = (pos / chr17_size) * 2 * np.pi
        radius = qual / 250 * 8 + 2

        if var_type == 'SNP':
            ax.plot(angle, radius, 'o', color='#3498db', markersize=qual / 40,
                    alpha=0.7, markeredgecolor='white', markeredgewidth=0.5)
        else:
            ax.plot(angle, radius, '^', color='#e74c3c', markersize=qual / 30,
                    alpha=0.8, markeredgecolor='white', markeredgewidth=0.5)

    # Chromosome backbone
    theta = np.linspace(0, 2 * np.pi, 360)
    ax.plot(theta, [1] * 360, color='#2c3e50', linewidth=3, alpha=0.3)

    # BRCA1 gene annotation
    brca1_pos = 43044295
    angle_brca1 = (brca1_pos / chr17_size) * 2 * np.pi
    ax.plot(angle_brca1, 8.5, 's', color='#f39c12', markersize=15)
    ax.text(angle_brca1, 9.5, 'BRCA1', ha='center', fontsize=10, fontweight='bold', color='#e67e22')

    ax.set_title(f'Circos Plot - {sample_name}\nChromosome 17',
                 fontsize=16, fontweight='bold', pad=30)
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498db',
               markersize=12, label=f'SNPs ({types.count("SNP")})'),
        Line2D([0], [0], marker='^', color='w', markerfacecolor='#e74c3c',
               markersize=12, label=f'INDELs ({types.count("INDEL")})'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#f39c12',
               markersize=12, label='BRCA1 Gene')
    ]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.3, 1.1))

    plt.tight_layout()
    output_file = f"{output_dir}/circos_plot.png"
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close()

    print(f"   Circos Plot: {output_file}")


def generate_plots(vcf_file, output_dir, sample_name):
    """
    Generate all plots: Manhattan, Heatmap, Circos

    Parameters:
    - vcf_file: Path to VCF file
    - output_dir: Directory to save plots
    - sample_name: Sample name for titles
    """
    print(f"\n{'=' * 55}")
    print("GENERATING PLOTS")
    print(f"{'=' * 55}\n")

    manhattan_plot(vcf_file, output_dir, sample_name)
    mutation_heatmap(vcf_file, output_dir, sample_name)
    circos_plot(vcf_file, output_dir, sample_name)

    print(f"\nAll plots generated in: {output_dir}")