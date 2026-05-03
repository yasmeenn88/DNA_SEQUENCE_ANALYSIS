#!/usr/bin/env python3
"""
ClinVar Validation Module
Converts genomic coordinates to HGVS format first, then searches ClinVar
"""

import requests
import json
import time


def convert_to_hgvs(chrom, pos, ref, alt):
    """
    Convert genomic coordinates to HGVS using NCBI E-utilities.
    Example: chr17:41234404 T>C → NM_007294.4(BRCA1):c.3113A>G
    """
    search_term = f"{chrom}:{pos}{ref}>{alt}"
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    try:
        # Step 1: Search ClinVar for the variant
        params = {
            'db': 'clinvar',
            'term': search_term,
            'retmax': 1,
            'retmode': 'json'
        }
        response = requests.get(url, params=params, timeout=10)
        data = response.json()
        id_list = data.get('esearchresult', {}).get('idlist', [])

        if not id_list:
            return None  # Variant not found

        # Step 2: Get summary to extract HGVS name
        variant_id = id_list[0]
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        fetch_params = {
            'db': 'clinvar',
            'id': variant_id,
            'retmode': 'json'
        }

        fetch_response = requests.get(fetch_url, params=fetch_params, timeout=10)
        fetch_data = fetch_response.json()
        result = fetch_data.get('result', {}).get(variant_id, {})

        # The title often contains the HGVS name
        title = result.get('title', '')

        # Try to find protein-level HGVS (starts with "NM_")
        if 'NM_' in title:
            return title

        # If not, return genomic HGVS
        return f"{chrom}:g.{pos}{ref}>{alt}"

    except Exception as e:
        return None


def search_clinvar_hgvs(hgvs_term):
    """
    Search ClinVar using HGVS format.
    Example: NM_007294.4(BRCA1):c.3113A>G
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        'db': 'clinvar',
        'term': hgvs_term,
        'retmax': 1,
        'retmode': 'json'
    }

    try:
        response = requests.get(url, params=params, timeout=10)
        data = response.json()
        id_list = data.get('esearchresult', {}).get('idlist', [])

        if not id_list:
            return {
                'hgvs_term': hgvs_term,
                'in_clinvar': False,
                'clinvar_id': None,
                'clinical_significance': 'Not in ClinVar'
            }

        # Get details
        variant_id = id_list[0]
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        fetch_params = {
            'db': 'clinvar',
            'id': variant_id,
            'retmode': 'json'
        }

        fetch_response = requests.get(fetch_url, params=fetch_params, timeout=10)
        fetch_data = fetch_response.json()
        result = fetch_data.get('result', {}).get(variant_id, {})

        return {
            'hgvs_term': hgvs_term,
            'in_clinvar': True,
            'clinvar_id': variant_id,
            'clinical_significance': result.get('clinical_significance', 'Not specified'),
            'title': result.get('title', '')
        }

    except Exception as e:
        return {
            'hgvs_term': hgvs_term,
            'in_clinvar': False,
            'error': str(e)
        }


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


def validate_with_clinvar(vcf_file, output_dir, max_variants=20):
    """
    Validate variants against ClinVar database using direct search.
    """
    print(f"\n{'='*55}")
    print("CLINVAR VALIDATION")
    print(f"{'='*55}")

    variants = parse_vcf(vcf_file)
    total = len(variants)
    to_check = min(total, max_variants)

    print(f"Total variants in VCF: {total}")
    print(f"Checking first {to_check} variants...\n")

    results = []
    in_clinvar_count = 0

    for i, var in enumerate(variants[:to_check]):
        genomic = f"{var['chrom']}:{var['pos']}{var['ref']}>{var['alt']}"
        print(f"  [{i+1}/{to_check}] {genomic}...", end=" ")

        # Search ClinVar directly
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            'db': 'clinvar',
            'term': f"{var['chrom']}[Chromosome] AND {var['pos']}[Position]",
            'retmax': 3,
            'retmode': 'json'
        }

        try:
            response = requests.get(url, params=params, timeout=10)
            data = response.json()
            id_list = data.get('esearchresult', {}).get('idlist', [])

            if id_list:
                in_clinvar_count += 1
                # Get clinical significance
                variant_id = id_list[0]
                fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                fetch_params = {'db': 'clinvar', 'id': variant_id, 'retmode': 'json'}
                fetch_response = requests.get(fetch_url, params=fetch_params, timeout=10)
                fetch_data = fetch_response.json()
                result = fetch_data.get('result', {}).get(variant_id, {})
                clinical_sig = result.get('clinical_significance', 'Unknown')

                results.append({
                    **var,
                    'in_clinvar': True,
                    'clinvar_id': variant_id,
                    'clinical_significance': clinical_sig,
                    'title': result.get('title', '')
                })
                print(f"FOUND - {clinical_sig}")
            else:
                results.append({
                    **var,
                    'in_clinvar': False,
                    'clinvar_id': None,
                    'clinical_significance': 'Not in ClinVar'
                })
                print("Not in ClinVar")

        except Exception as e:
            results.append({
                **var,
                'in_clinvar': False,
                'error': str(e)
            })
            print(f"Error: {e}")

        time.sleep(0.3)

    # Save results
    import json
    output_file = f"{output_dir}/clinvar_validation.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    # Summary
    print(f"\n{'='*55}")
    print("CLINVAR VALIDATION SUMMARY")
    print(f"{'='*55}")
    print(f"   Total checked: {to_check}")
    print(f"   Found in ClinVar: {in_clinvar_count}/{to_check}")

    pathogenic = sum(1 for r in results if 'pathogenic' in str(r.get('clinical_significance', '')).lower())
    if pathogenic > 0:
        print(f"   Pathogenic/Likely Pathogenic: {pathogenic}")

    print(f"   Results saved: {output_file}")
    print(f"{'='*55}\n")

    return results