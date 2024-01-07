#!/usr/bin/env python3

import argparse, sys, re, os
import pysam

__version__ = '1.0.0'

def checkfile(file_path):
    """Check if a file exists at the given path."""
    if not os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    return file_path

def fileexists(file_path):
    """Check if a file exists at the given path."""
    if os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The outfile {file_path} exists!")
    return file_path

def natural_sort_key(s):
    """A sorting key for natural string order."""
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

def main():

    parser = argparse.ArgumentParser(description="Merge small variants and gene-level SVs and apply filtering criteria")
    parser.add_argument("smallvariantvcffile",  type=checkfile, help="Small variant VCF file")
    parser.add_argument("svvcffile", type=checkfile, help="SV VCF file")
    parser.add_argument("-b", "--bedfile", type=checkfile, required=True, help="Bed file to filter VCFs")
    parser.add_argument("-L", "--maxsvlength", type=int, default=500, help="Maximum gene-level SV length")
    parser.add_argument("-m", "--minreads", type=int, default=3, help="Minimum alt-supporting reads to pass LowReads filter")
    parser.add_argument("-a", "--minvaf", type=float, default=5.0, help="Minimum VAF to pass MinVAF filter")
    parser.add_argument("-o", "--outfile", type=fileexists, help="Outfile")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s: ' + __version__)

    args = parser.parse_args()

    outfile = args.outfile
    if outfile is None:
        outfile = sys.stdout

    # Open VCF file
    vcf = pysam.VariantFile(args.smallvariantvcffile, "r")
    svVcf = pysam.VariantFile(args.svvcffile, "r")

    header = vcf.header.copy()
    vcf_out = pysam.VariantFile(outfile, "w", header=header)
    vcf_out.header.merge(svVcf.header)

    if 'LowReads' not in vcf_out.header.filters.keys():
        vcf_out.header.filters.add("LowReads",None,None,f'Fails minimum alt-supporting reads (${args.minreads})')

    if 'MinVAF' not in vcf_out.header.filters.keys():
        vcf_out.header.filters.add("MinVAF",None,None,f'Fails minimum VAF filter (${args.minvaf})')

    regions = pysam.TabixFile(args.bedfile)

    variants = []

    for row in regions.fetch(parser=pysam.asTuple()):
        reg = f'{row[0]}:{row[1]}-{row[2]}'
        
        for record in vcf.fetch(region=reg):

            new_record = vcf_out.new_record()

            new_record.chrom = record.chrom
            new_record.pos = record.pos
            new_record.stop = record.stop
            new_record.id = record.id
            new_record.alleles = record.alleles

            for k in record.filter.keys():
                new_record.filter.add(k)

            for k in record.info.keys():
                new_record.info[k] = record.info.get(k)

            for k in record.format.keys():
                new_record.samples[0][k] = record.samples[0][k]

            if record.samples[0]['AD'][1] < int(args.minreads):
                new_record.filter.add("LowReads")

            if record.samples[0]['AF'][0] * 100 < float(args.minvaf):
                new_record.filter.add("MinVAF")

            variants = variants + [new_record]

        for record in svVcf.fetch(region=reg):

            # skip large SV records and  
            if record.info['SVTYPE'] == 'BND' or ('SVLEN' in record.info.keys() and abs(record.info['SVLEN'][0]) > int(args.maxsvlength)):
                continue

            new_record = vcf_out.new_record()

            new_record.chrom = record.chrom
            new_record.pos = record.pos
            new_record.stop = record.stop
            new_record.id = record.id
            new_record.alleles = record.alleles

            for k in record.filter.keys():
                new_record.filter.add(k)

            for k in record.info.keys():
                new_record.info[k] = record.info.get(k)

            for k in record.format.keys():
                new_record.samples[0][k] = record.samples[0][k]


            PR = (0,0)
            SR = (0,0)
            if 'PR' in record.format.keys():
                PR = record.samples[0]['PR']

            if 'SR' in record.format.keys():
                SR = record.samples[0]['SR']

            DP = PR[0] + PR[1] + SR[0] + SR[1]
            AD = (PR[0] + SR[0], PR[1] + SR[1])
            AF = round((PR[1]+SR[1]) / DP,4)

            new_record.samples[0]['DP'] = DP
            new_record.samples[0]['AD'] = AD
            new_record.samples[0]['AF'] = AF

            if new_record.samples[0]['AD'][1] < int(args.minreads):
                new_record.filter.add("LowReads")

            if new_record.samples[0]['AF'][0] * 100 < float(args.minvaf):
                new_record.filter.add("MinVAF")

            variants = variants + [new_record]

    for record in variants:
        vcf_out.write(record)

    vcf_out.close()

if __name__ == "__main__":
    main()
