#!/usr/bin/env python3
import argparse, re, sys, itertools, os
import pandas as pd
import pyranges as pr
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
        raise argparse.ArgumentTypeError(f"The ouytfile {file_path} exists!")
    return file_path

def ratio2abundance(cn, nl, l2r):
    t = (nl * (2**l2r - 1)) / (cn - nl)
    return t

def process_ichor(segsfile,gender):
    records = pd.DataFrame(columns=["Chromosome", "Start", "End", "ID", "REF", "ALT","CN","QUAL", "FILTER","INFO","FORMAT"])
    with open(segsfile, "r") as CNV:
        next(CNV)  # Skip header
        for line in CNV:
            line = line.strip()
            F = line.split("\t")

            if gender == "male" and F[1] == "chrY" and float(F[5]) < -0.38:
                F[7] = "HETD"
                F[6] = str(int(F[6]) - 1)

            if "NEUT" in F[7]:
                continue
            
            svtype = "DEL" if float(F[5]) <= 0 else "DUP"
            pos = int(F[2]) - 1
            refnt = 'N'
            svlen = abs(int(F[3]) - int(F[2]) + 1)

            nl = 2
            if bool(re.search(r"[X|Y]", F[1])) and gender == "male":
                nl = 1
            
            if float(F[5]) < 0:
                F[6] = str(nl - 1)
            else:
                F[6] = str(nl + 1)

            abund = round(ratio2abundance(int(F[6]), nl, float(F[5])), 2)
            
            filter = 'PASS'
            
            if svtype == "DEL":
                svlen = -svlen
            
            info_dict = dict(zip(['SVTYPE','SVLEN','set'],['CNV',svlen,'IchorCNA']))
            format_dict = dict(zip(['GT','BN','CN','CF'],[(0,1),int(F[4]),int(F[6]),abund]))
            
            records.loc[len(records)] = [F[1], str(pos), F[3], f"ICHOR:{F[1]}_{F[2]}_{F[3]}", refnt,
                                f"<{svtype}>",int(F[6]), ".", filter, info_dict, format_dict]
            
    return pr.PyRanges(records)

def vcf_to_pyranges(vcf_file):
    vcf = pysam.VariantFile(vcf_file)

    records = pd.DataFrame(columns=["Chromosome", "Start", "End", "ID", "REF", "ALT","CN","QUAL","FILTER","INFO","FORMAT"])
    for record in vcf.fetch():
        chrom = record.contig
        start = record.pos
        end = record.stop
        record_id = record.id
        ref = record.ref
        alt = ','.join(record.alts)
        qual = record.qual
        filter = ';'.join(record.filter.keys())
        
        info_dict = dict(record.info)
        info_dict['set'] = 'dragen'
        format_dict = {}
        for field in vcf.header.formats.keys():
            if field in record.samples[0].keys():
                format_dict[field] = record.samples[0][field]
            else:
                format_dict[field] = None

        row_data = [chrom, start, end, record_id, ref, alt, format_dict['CN'], qual, filter, info_dict,format_dict]

        records.loc[len(records)] = row_data

    return pr.PyRanges(records)


def merge_records(df,overlap=0.9):

    if df is None or df.empty:
        return pd.DataFrame(columns=df.columns) 
    
    if 'Cluster' in df.columns:
        df = df.drop(columns='Cluster')

    # unmerged records
    unmergedDf = pd.DataFrame(columns=df.columns)
    # merged records
    mergedDf = pd.DataFrame(columns=df.columns)
    
    # if there is only one record then done merge it.
    if df.shape[0] == 1:
        unmergedDf = df.copy()

    else:
        # records to merge
        to_merge = pd.DataFrame(columns=df.columns)

        # get each record and determine overlap with remaining records
        for index, item in df.iterrows():            
            row = pr.PyRanges(df[df.index == index])
            cluster = pr.PyRanges(df[df.index != index]).merge()
            ov1 = row.coverage(cluster).df.to_dict()['FractionOverlaps'][0]
            ov2 = cluster.coverage(row).df.to_dict()['FractionOverlaps'][0]
            
            # if either the record or the cluster has an overlap greater than the threshold then merge them
            if ov1 > overlap or ov2 > overlap:
                to_merge = pd.concat([to_merge,row.df],axis=0)
            else: # otherwise add to unmerged records
                unmergedDf = pd.concat([unmergedDf,row.df],axis=0)

        # if there are records to merge then merge them
        if to_merge.shape[0] > 0:
            merged = pr.PyRanges(to_merge).merge().df.loc[0].to_dict()
            merged['ID'] = ';'.join(to_merge['ID'].to_list())
            merged['REF'] = 'N'
            merged['ALT'] = to_merge['ALT'].tolist()[0]
            merged['QUAL'] = pd.to_numeric(to_merge['QUAL'],errors='coerce').max()
            if to_merge[to_merge['FILTER'].str.contains('PASS')].shape[0] > 0:
                merged['FILTER'] = 'PASS'
            else:
                merged['FILTER'] = ';'.join(to_merge['FILTER'].to_list())
            
            info_dict = {'SVTYPE':'CNV','SVLEN':merged['End'] - merged['Start'] + 1}
            
            if merged['ALT'] == '<DEL>' or merged['ALT'] == '<DEL>,<DUP>' or merged['ALT'] == '<DUP>,<DEL>':
                info_dict['REFLEN'] = info_dict['SVLEN']
            else :
                info_dict['REFLEN'] = -info_dict['SVLEN']

            if merged['ALT'] == '<DEL>':
                info_dict['SVLEN'] = -info_dict['SVLEN']

            format_dict = {'GT':(0,1),'CN':int(to_merge['CN'].tolist()[0])}
            # get abundance estimate. This is the average of the abundance from ICHOR and 2 * the MAF from Dragen.
            abundance = []
            callers = []
            for index, item in to_merge.iterrows():
                if 'MAF' in item['FORMAT'].keys():
                    if item['FORMAT']['MAF'] is not None:
                        abundance = abundance + [item['FORMAT']['MAF'] * 2]
                elif 'CF' in item['FORMAT'].keys():
                    abundance = abundance + [item['FORMAT']['CF']]
                callers = callers + [item['INFO']['set']]

            info_dict['set'] = ','.join(callers)

            if len(abundance) > 0:
                format_dict['CF'] = round(sum(abundance) / len(abundance),1) * 100
            else:
                format_dict['CF'] = None

            mergedDf = pd.DataFrame([merged])
            mergedDf['INFO'] = [info_dict]
            mergedDf['FORMAT'] = [format_dict]

    # rescale Dragen's MAF to a AF for the unmerged records
    for rec in unmergedDf.to_dict(orient='records'):
        if 'CF' not in rec['FORMAT'].keys() and 'MAF' in rec['FORMAT'].keys() and rec['FORMAT']['MAF'] is not None:
            rec['FORMAT']['CF'] = rec['FORMAT']['MAF'] * 2 * 100
        elif 'CF' not in rec['FORMAT'].keys():
            rec['FORMAT']['CF'] = None

        mergedDf = pd.concat([mergedDf,pd.DataFrame([rec])],axis=0)

    return mergedDf


def main():
    parser = argparse.ArgumentParser(description="Merge VCF records with overlap")
    parser.add_argument("ichorsegfile", type=checkfile, help="Path to the first input VCF file")
    parser.add_argument("vcffile", type=checkfile, help="Path to the second input VCF file")
    parser.add_argument("-c","--cytobandbed", type=checkfile, help="Bed file with cytoband intervals")
    parser.add_argument("-g", "--gender", required=True, choices=["male", "female"], help="Gender ('male' or 'female')")
    parser.add_argument("-s", "--minsize", type=int, default=5000000, help="Minimum CNA size")
    parser.add_argument("-f", "--minabund", type=float, default=10.0, help="Minimum CNA cell fraction to pass filter")
    parser.add_argument("-m", "--discardabund", type=float, default=0.05, help="Discard abundance")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    # process Ichor input
    ichorPr = process_ichor(args.ichorsegfile,args.gender)
    vcfPr = vcf_to_pyranges(args.vcffile)

    mergedCnvDf = pr.PyRanges(pd.concat([ichorPr.df,vcfPr.df],axis=0)).cluster(by=['ALT','CN'],slack=-1).df
    mergedCnvDf = mergedCnvDf.groupby('Cluster').apply(merge_records)

    # make new VCF and set headers
    inVcf = pysam.VariantFile(args.vcffile)

    # add info tags
    if 'set' not in inVcf.header.info.keys():
        inVcf.header.info.add("set",1,'String','Callers that identified this variant')

    if 'CYTOBANDS' not in inVcf.header.info.keys():
        inVcf.header.info.add("CYTOBANDS",1,'String','List of cytobands in order that overlap the variant')

    if 'BN' not in inVcf.header.formats.keys():
        inVcf.header.formats.add("BN", 1, 'Integer', 'IchorCNA bins spanning this variant')

    if 'CF' not in inVcf.header.formats.keys():
        inVcf.header.formats.add("CF", 1, 'Float', 'Estimated tumor cell fraction of variant (%)')

    # add filters
    if 'MinCNVSize' not in inVcf.header.filters.keys():
        inVcf.header.filters.add("MinCNVSize",None,None,'CNV does not meet minimum size of '+str(args.minsize)+' bp')

    # add filters
    if 'MinCNVAbundance' not in inVcf.header.filters.keys():
        inVcf.header.filters.add("MinCNVAbundance",None,None,'CNV does not meet minimum tumor cell fraction of '+str(args.minabund))

    # Create a new VCF writer with added header lines
    outVcfFile = sys.stdout
    outVcf = pysam.VariantFile(outVcfFile, 'w', header=inVcf.header)

    cytobandsPr = None
    if args.cytobandbed is not None:
        cytobandsPr = pr.PyRanges(pd.read_csv(args.cytobandbed, header=None, names="Chromosome Start End Band Color".split(), sep="\t"))

    for item in mergedCnvDf.to_dict(orient='records'):
        nrec = outVcf.new_record()
        nrec.chrom = item['Chromosome']
        nrec.pos = item['Start']
        nrec.stop = item['End']
        nrec.id = item['ID']
        nrec.filter.clear()

        # add existing filters
        if 'PASS' in item['FILTER']:
            nrec.filter.add('PASS')
        else:
            filters = item['FILTER'].split(';')
            if len(set(['cnvLength','cnvQual','binCount']) & set(filters))>0:
                continue
            for v in filters:
                nrec.filter.add(v)

        # add new filters        
        if item['End'] - item['Start'] < int(args.minsize):
            if 'PASS' not in nrec.filter.keys():
                nrec.filter.add('MinCNVSize')
            else:
                nrec.filter.clear()
                nrec.filter.add('MinCNVSize')

        if item['FORMAT']['CF'] is None or item['FORMAT']['CF'] < float(args.minabund):
            if 'PASS' not in nrec.filter.keys():
                nrec.filter.add('MinCNVAbundance')
            else:
                nrec.filter.clear()
                nrec.filter.add('MinCNVAbundance')

        nrec.ref = item['REF']
        nrec.alts = (item['ALT'],)

        for k in inVcf.header.info.keys():
            if k == 'END':
                continue
            if k in item['INFO'].keys():
                nrec.info[k] = item['INFO'].get(k)
            else:
                nrec.info[k] = None

        # add cytobands, if present
        if cytobandsPr is not None:
            bands = '&'.join(cytobandsPr.intersect(pr.PyRanges(chromosomes = str(item['Chromosome']),starts = [item['Start']],ends = [item['End']])).df['Band'].tolist())
            nrec.info['CYTOBANDS'] = bands

        for k in inVcf.header.formats.keys():
            if k in item['FORMAT'].keys():
                nrec.samples[0][k] = item['FORMAT'][k]
            else:
                nrec.samples[0][k] = tuple(itertools.repeat(None, inVcf.header.formats.get(k).number))

        outVcf.write(nrec)

    outVcf.close()

if __name__ == "__main__":
    main()
