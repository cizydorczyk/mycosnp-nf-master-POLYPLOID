#!/usr/bin/env python
"""
Filter genotypes of a VCF using Genotype Quality Score (GQ), Depth (DP),
Allelic Depth (AD) and Allelic ratio via binomial test.
Handles both variant sites (SNPs) and invariant reference sites.
"""


from __future__ import print_function
from __future__ import division


import sys
import re
import argparse
import vcfTools


# take input
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('infile', help='VCF filename', type=str)
parser.add_argument('--min_GQ', help='minimum GQ/RGQ to be kept', type=int)
parser.add_argument('--keep_GQ_0_refs', help='keep ref calls with GQ/RGQ 0 despite min_GQ', action='store_true')
parser.add_argument('--min_percent_alt_in_AD', help='min percent alt in AD to be kept for variants', type=float)
parser.add_argument('--min_total_DP', help='min total DP to be kept', type=float)
parser.add_argument('--het_binomial_p', help='filter hets that fail binomial test with given p-value', type=float)
parser.add_argument('--keep_all_ref', help='don\'t filter reference bases', action='store_true')
args = parser.parse_args()



# set default parameters
infile = args.infile
keep_all_ref = args.keep_all_ref
keep_GQ_0_refs = args.keep_GQ_0_refs
min_GQ = int(0)
min_AD = float(0)
min_tot_DP = int(0)
het_binomial_p = False


if args.min_GQ:
    min_GQ = args.min_GQ
if args.min_percent_alt_in_AD:
    min_AD = args.min_percent_alt_in_AD
if args.min_total_DP:
    min_tot_DP = args.min_total_DP
if args.het_binomial_p:
    het_binomial_p = args.het_binomial_p


# initialize analysis
comment_pattern = re.compile(r"^#")


genome_list = []


header = vcfTools.VcfHeader(infile)
caller = header.get_caller()
samples = header.get_samples()


# initialize summary statistics
filtered_GQ = {}
for sample in samples:
    filtered_GQ[sample] = 0


filtered_AD = {}
for sample in samples:
    filtered_AD[sample] = 0


filtered_DP = {}
for sample in samples:
    filtered_DP[sample] = 0


filtered_Bi = {}
for sample in samples:
    filtered_Bi[sample] = 0


filtered = {}
for sample in samples:
    filtered[sample] = 0



# Helper function to check if genotype is reference
def is_ref_genotype(gt):
    return gt in ['0', '0/0', '0|0', '0/0/0', '0|0|0']


# Helper function to check if genotype is missing
def is_missing_genotype(gt):
    return gt in ['.', './.', '.|.', '././.', '.|.|.']


# Helper function to check if a site is a variant site (has ALT alleles)
def is_variant_site(record):
    """Check if site is a variant (SNP/indel) or invariant (reference-only)"""
    alt_field = record.get_alt_field()
    # Invariant sites typically have '.' as ALT or no proper ALT allele
    return alt_field != '.' and alt_field != ''


# Helper function to extract RGQ from genotype field for reference sites
def get_RGQ_from_genotype(genotype_field, format_string):
    """Extract RGQ value from genotype field for invariant sites"""
    format_fields = format_string.split(':')
    genotype_values = genotype_field.split(':')
    
    if 'RGQ' in format_fields:
        rgq_index = format_fields.index('RGQ')
        if rgq_index < len(genotype_values):
            rgq_value = genotype_values[rgq_index]
            return rgq_value if rgq_value != '.' else None
    return None


# Helper function to extract DP from genotype field for reference sites
def get_DP_from_genotype(genotype_field, format_string):
    """Extract DP value from genotype field"""
    format_fields = format_string.split(':')
    genotype_values = genotype_field.split(':')
    
    if 'DP' in format_fields:
        dp_index = format_fields.index('DP')
        if dp_index < len(genotype_values):
            dp_value = genotype_values[dp_index]
            return dp_value if dp_value != '.' else None
    return None


# Helper function to extract AD from genotype field for reference sites
def get_AD_from_genotype(genotype_field, format_string):
    """Extract AD value from genotype field for reference sites"""
    format_fields = format_string.split(':')
    genotype_values = genotype_field.split(':')
    
    if 'AD' in format_fields:
        ad_index = format_fields.index('AD')
        if ad_index < len(genotype_values):
            ad_value = genotype_values[ad_index]
            return ad_value if ad_value != '.' else None
    return None


# Helper function to filter reference site genotypes
def filter_reference_genotype(genotype_field, format_string, min_rgq, min_per_ad, min_tot_dp, keep_rgq_0):
    """
    Filter reference site genotypes based on RGQ, AD, and DP.
    Mimics the logic used for variant sites but adapted for invariant site FORMAT fields.
    Returns: (filtered_genotype_field, RGQ_flag, AD_flag, DP_flag)
    """
    split_genotype = genotype_field.split(':')
    gt = split_genotype[0]
    
    RGQ_flag = False
    AD_flag = False
    DP_flag = False
    
    # Extract RGQ, AD, and DP values
    rgq = get_RGQ_from_genotype(genotype_field, format_string)
    ad = get_AD_from_genotype(genotype_field, format_string)
    dp = get_DP_from_genotype(genotype_field, format_string)
    
    # Convert to numeric values for comparison
    rgq_numeric = None
    if rgq is not None:
        try:
            rgq_numeric = int(rgq)
        except ValueError:
            try:
                rgq_numeric = float(rgq)
            except ValueError:
                pass
    
    ad_numeric = None
    if ad is not None:
        try:
            ad_numeric = int(ad)
        except ValueError:
            try:
                ad_numeric = float(ad)
            except ValueError:
                pass
    
    dp_numeric = None
    if dp is not None:
        try:
            dp_numeric = int(dp)
        except ValueError:
            try:
                dp_numeric = float(dp)
            except ValueError:
                pass
    
    # Apply filters
    should_filter = False
    
    # Handle RGQ=0 special case if keep_GQ_0_refs is set
    if keep_rgq_0 and rgq_numeric is not None and rgq_numeric == 0:
        # Only check AD and DP for RGQ=0 refs when keep_GQ_0_refs is set
        # Check AD percent (allelic balance for reference)
        if ad_numeric is not None and dp_numeric is not None and dp_numeric > 0:
            ad_percent = (ad_numeric / dp_numeric) * 100
            if ad_percent < min_per_ad:
                should_filter = True
                AD_flag = True
        
        # Check DP
        if dp_numeric is not None and dp_numeric < min_tot_dp:
            should_filter = True
            DP_flag = True
    else:
        # Normal filtering logic - apply all filters
        # Check RGQ
        if rgq_numeric is not None and rgq_numeric < min_rgq:
            should_filter = True
            RGQ_flag = True
        
        # Check AD percent (allelic balance for reference)
        if ad_numeric is not None and dp_numeric is not None and dp_numeric > 0:
            ad_percent = (ad_numeric / dp_numeric) * 100
            if ad_percent < min_per_ad:
                should_filter = True
                AD_flag = True
        
        # Check DP
        if dp_numeric is not None and dp_numeric < min_tot_dp:
            should_filter = True
            DP_flag = True
    
    # If filters fail, set genotype to missing
    if should_filter:
        split_genotype[0] = './.'
        return ':'.join(split_genotype), RGQ_flag, AD_flag, DP_flag
    
    return genotype_field, RGQ_flag, AD_flag, DP_flag



# process vcf and apply genotype filters
with open(infile, 'r') as vcf_file:
    for vcf_line in vcf_file:
        # output header
        if (re.search(comment_pattern, vcf_line)):
            print(vcf_line, end="")


        else:
            # output site level information
            record = vcfTools.VcfRecord(vcf_line)


            # Skip records with deletion alleles (*)
            #ref = record.get_ref()
            #alt_field = record.get_alt_field()
            
            #if ref == '*' or '*' in alt_field:
            #   continue  # Skip this entire record - don't output it


            new_vcf_line = "\t".join([str(record.get_chrom()),
                                     str(record.get_pos()),
                                     str(record.get_id()),
                                     str(record.get_ref()),
                                     str(record.get_alt_field()),
                                     str(record.get_qual()),
                                     str(record.get_filter()),
                                     str(record.get_info()),
                                     str(record.get_format())])
            print(new_vcf_line, end="")
            
            
            # Determine if this is a variant site or reference site
            is_variant = is_variant_site(record)


            # filter genotypes
            for sample in samples:
                genotype_field = record.get_genotypes_field(header.get_sample_index(sample))
                split_genotype = genotype_field.split(':')
                init_gt = split_genotype[0]
                
                # Process variant sites using existing vcfTools logic
                if is_variant:
                    gq = record.get_GQ(split_genotype[0], index=header.get_sample_index(sample))
                    
                    if is_ref_genotype(split_genotype[0]) and keep_all_ref:
                        # get_genotype_REF() allows testing of ref AD support against 80% threshold
                        split_genotype[0], GQ_flag, AD_flag, DP_flag, Bi_flag = record.get_genotype_REF(
                            index=header.get_sample_index(sample), return_flags=True)
                    elif is_ref_genotype(split_genotype[0]) and (gq == '0') and keep_GQ_0_refs:
                        split_genotype[0], GQ_flag, AD_flag, DP_flag, Bi_flag = record.get_genotype_REF(
                            index=header.get_sample_index(sample), min_per_ad=80.0, min_tot_dp=min_tot_DP, return_flags=True)
                    elif is_ref_genotype(split_genotype[0]):
                        split_genotype[0], GQ_flag, AD_flag, DP_flag, Bi_flag = record.get_genotype_REF(
                            index=header.get_sample_index(sample), min_gq=min_GQ, min_per_ad=80.0, min_tot_dp=min_tot_DP, return_flags=True)
                    else:
                        # Variant, so test with get_genotype2()
                        split_genotype[0], GQ_flag, AD_flag, DP_flag, Bi_flag = record.get_genotype2(
                            index=header.get_sample_index(sample), min_gq=min_GQ, min_per_ad=min_AD,
                            min_tot_dp=min_tot_DP, het_binom_p=het_binomial_p, return_flags=True)
                    
                    # aggregate summary statistics
                    if not is_missing_genotype(init_gt):
                        if GQ_flag:
                            filtered_GQ[sample] += 1
                        if AD_flag:
                            filtered_AD[sample] += 1
                        if DP_flag:
                            filtered_DP[sample] += 1
                        if Bi_flag:
                            filtered_Bi[sample] += 1
                        if GQ_flag or AD_flag or DP_flag or Bi_flag:
                            filtered[sample] += 1
                    
                    print("\t".join(["", str(":".join(split_genotype))]), end="")
                
                # Process invariant (reference only) sites using new reference site logic
                else:
                    format_string = record.get_format()
                    
                    # Apply reference site filtering with same logic as variants
                    if keep_all_ref:
                        # Don't filter at all
                        filtered_genotype_field = genotype_field
                        GQ_flag = False
                        AD_flag = False
                        DP_flag = False
                        Bi_flag = False
                    else:
                        # Apply same filters as variant sites
                        # Use filter_reference_genotype() because this is an INVARIANT site, not just a ref genotype
                        filtered_genotype_field, GQ_flag, AD_flag, DP_flag = filter_reference_genotype(
                            genotype_field, format_string, min_rgq=min_GQ, min_per_ad=80.0, # Hard-coded 80% support for ref base required!
                            min_tot_dp=min_tot_DP, keep_rgq_0=keep_GQ_0_refs)
                        Bi_flag = False  # No binomial test for reference sites
                    
                    # aggregate summary statistics for reference sites
                    if not is_missing_genotype(init_gt):
                        if GQ_flag:
                            filtered_GQ[sample] += 1
                        if AD_flag:
                            filtered_AD[sample] += 1
                        if DP_flag:
                            filtered_DP[sample] += 1
                        if GQ_flag or AD_flag or DP_flag:
                            filtered[sample] += 1
                    
                    print("\t".join(["", filtered_genotype_field]), end="")


            # output end of VCFs
            print("")


# output summary statistics
print("Sample", file=sys.stderr, end="\t")
print("\t".join(samples), file=sys.stderr)
print("filtered_GQ", file=sys.stderr, end="")
for sample in samples:
    print("".join(["\t", str(filtered_GQ[sample])]), end="", file=sys.stderr)
print("\nfiltered_AD", file=sys.stderr, end="")
for sample in samples:
    print("".join(["\t", str(filtered_AD[sample])]), end="", file=sys.stderr)
print("\nfiltered_DP", file=sys.stderr, end="")
for sample in samples:
    print("".join(["\t", str(filtered_DP[sample])]), end="", file=sys.stderr)
print("\nfiltered_Bi", file=sys.stderr, end="")
for sample in samples:
    print("".join(["\t", str(filtered_Bi[sample])]), end="", file=sys.stderr)
print("\nfiltered_total", file=sys.stderr, end="")
for sample in samples:
    print("".join(["\t", str(filtered[sample])]), end="", file=sys.stderr)
print("", file=sys.stderr)
