#!/usr/bin/env python3
import argparse, sys, os, csv
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

# write a function to parse the VEP CSQ VCF header
def vepHeader(header):
    fields = header.info.get('CSQ').record.get('Description').split(':')[1].replace('"',"").strip().split('|')
    return dict(zip(fields,range(0,len(fields))))

def parseVepCsq(csq, vep, field=None):
    if not isinstance(csq, tuple):
        raise TypeError("csq should be a string")
    if not isinstance(vep, list):
        raise TypeError("vep should be a list")
    
    out = []
    for val in csq:
        fields = val.split("|")
        
        if len(fields) < len(vep):
            raise ValueError("VEP header has more fields than the CSQ record")
        
        csq_dict = dict(zip(vep, fields))
        out.append(csq_dict)

    if field:
        out = list(set([d[field] for d in out if field in d and d[field]!='']))

    return out

def main():
    
    parser = argparse.ArgumentParser(description="Filter VCF records based on criteria")
    parser.add_argument("vcffile", type=checkfile, help="Input VCF file")
    parser.add_argument("-g", "--geneList", type=checkfile, required=True, help="Gene list file")
    parser.add_argument("-l", "--minSvLength", type=int, default=5000, help="Minimum SV length")
    parser.add_argument("-L", "--maxSvLength", type=int, default=10000000, help="Maximum SV length")
    parser.add_argument("-m", "--minSvReads", type=int, default=2, help="Minimum SR and PR alt-supporting reads")
    parser.add_argument("-a", "--minSvAbundance", type=float, default=5.0, help="Minimum SV abundance in percent")
    parser.add_argument("-o", "--outfile", type=fileexists, help="Outfile")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s: ' + __version__)

    nonSynon = {"splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_ablation","transcript_amplification","inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant"}

    args = parser.parse_args()

    outfile = args.outfile
    if outfile is None:
        outfile = sys.stdout

    # Read the gene list and store unique genes in knownGenes
    knownGenes = set()
    csvfile = open(args.geneList, "r")
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        knownGenes.add(row[0])
        knownGenes.add(row[2])

    passingRecords = set()
    knownSvGenes = {}

    # Open VCF file
    vcf_in = pysam.VariantFile(args.vcffile, "r")

    vep = vepHeader(vcf_in.header)

    for record in vcf_in:
        # this gets genes that overlap SVs from both the standard VEP annotation to get upstream/downstream events and custom BED overlap annotations
        # so we dont miss anything. It also annotates variants with gene clusters, like IGH, IGL, TRA, etc.  
        genes = list(set(parseVepCsq(record.info.get("CSQ"),list(vep.keys()),'SYMBOL') + parseVepCsq(record.info.get("CSQ"),list(vep.keys()),'KnownSvGenes')))
        svtype = record.info.get("SVTYPE")

        # skip dels/dups that are small or very large--these are captured via other means
        if svtype in ['DEL','DUP'] and record.info.get("SVLEN") is not None:
            if abs(record.info.get("SVLEN")[0]) < args.minSvLength or abs(record.info.get("SVLEN")[0]) > args.maxSvLength:
                continue

        consequences = '&'.join(parseVepCsq(record.info.get("CSQ"),list(vep.keys()),'Consequence')).split('&')
    
        PR = (0,0)
        SR = (0,0)
        if 'PR' in record.format.keys():
            PR = record.samples[0]['PR']

        if 'SR' in record.format.keys():
            SR = record.samples[0]['SR']

        # some records have no read support. Skip these.
        if PR[1] == 0 and SR[1] == 0:
            continue

        # csq has a known gene
        recordVariant = False
        geneHits = set(genes) & knownGenes

        # get records that overlap a known gene
        if len(geneHits) > 0 and svtype=='BND':
            recordVariant = True

        if len(geneHits) > 0 and svtype in ['INS','DEL','DUP'] and len(set(consequences) & nonSynon)>0:
            recordVariant = True
            
        # also get BND records where the mate has been recorded (this is the main purpose of this script) 
        if svtype=='BND' and record.info.get("MATEID")[0] in passingRecords:
            recordVariant = True

        # finally get passing records with a functional consequence if there is a contig and there are at least minReads alternate-supporting SR and PR reads
        if len(record.filter)==0 and len(set(consequences) & nonSynon)>0 and record.info.get("CONTIG") is not None and PR[1] >= args.minSvReads and SR[1] >= args.minSvReads:
            recordVariant = True

        if recordVariant:
            passingRecords.add(record.id)
            if len(geneHits) > 0:
                knownSvGenes[record.id] = list(geneHits)

    vcf_in.close()

    # Re-open the VCF and print header and records that are in passingRecords
    vcf_in = pysam.VariantFile(args.vcffile, "r")

    if 'KnownSvGenes' not in vcf_in.header.info.keys():
        vcf_in.header.info.add("KnownSvGenes", '.', 'String', 'List of recurrent SV genes that overlap this variant')

    if 'Imprecise' not in vcf_in.header.filters.keys():
        vcf_in.header.filters.add("Imprecise",None,None,'No contig found so breakend are imprecise')

    if 'MinSvReads' not in vcf_in.header.filters.keys():
        vcf_in.header.filters.add("MinSvReads",None,None,f'Fails minimum SR or PR alt-supporting reads (${args.minSvReads})')

    if 'MinSvAbundance' not in vcf_in.header.filters.keys():
        vcf_in.header.filters.add("MinSvAbundance",None,None,f'Fails minimum SV abundance (${args.minSvAbundance})')

    header = vcf_in.header.copy()

    vcf_out = pysam.VariantFile(outfile, "w", header=header)

    for record in vcf_in:
        printRecord = False

        if record.id in passingRecords:
            printRecord = True
        elif record.info.get("SVTYPE") == 'BND' and record.info.get("MATEID")[0] in passingRecords:
            printRecord = True

        if printRecord:

            PR = (0,0)
            SR = (0,0)
            if 'PR' in record.format.keys():
                PR = record.samples[0]['PR']

            if 'SR' in record.format.keys():
                SR = record.samples[0]['SR']

            # apply filters
            if record.info.get("CONTIG") is None:
                record.filter.add('Imprecise')

            if PR[1] < int(args.minSvReads) or SR[1] < int(args.minSvReads):
                record.filter.add('MinSvReads')

            if (PR[1] + PR[1]) / (PR[1] + SR[1] + PR[0] + PR[1] + SR[0] + SR[1]) * 100.0 < float(args.minSvAbundance):
                record.filter.add('MinSvAbundance')

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

            if new_record.id in knownSvGenes.keys():
                new_record.info['KnownSvGenes'] = ','.join(knownSvGenes[new_record.id])

            for k in record.format.keys():
                new_record.samples[0][k] = record.samples[0][k]
            
            vcf_out.write(new_record)

    vcf_in.close()
    vcf_out.close()

if __name__ == "__main__":
    main()
