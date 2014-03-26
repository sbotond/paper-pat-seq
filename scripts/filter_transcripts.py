"""
Parse out the transcriptome using a reference genome and a GFF file. Filter transcripts 
by number of peaks.
"""
import HTSeq
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
from Bio import Alphabet
import argparse
import re
import pickle
import operator
from collections import defaultdict

chr_RE = re.compile("chr[IVX]+")
chr_gff_RE = re.compile("\\[chromosome=([IVX]+)\\]")
mito_gff_RE = re.compile("\\[location=mitochondrion\\]")

# Genomic data from
# http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/

roman = [ "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI" ]
int2roman = dict(zip(range(1, 16+1), roman))
roman2int = dict(zip(roman, range(1, 16+1)))

def parse_args():
    argparser = argparse.ArgumentParser(description="Filter transcripts on the basis GFF files and poly(A) peaks.")
    argparser.add_argument("-r", metavar="ref", type=str, help="FASTA file with the reference genome", required=True)
    argparser.add_argument("-t", metavar="tr", type=str, help="GFF with the coding transcriptome", required=True)
    argparser.add_argument("-p", metavar="peaks", type=str, help="GFF with the peaks", required=True)
    argparser.add_argument("-ot", metavar="out_fasta", type=str, help="Output FASTA file with the full transcriptome", \
                           required=True)
    argparser.add_argument("-o", metavar="out_filtered", type=str, help="Output file with the filtered transcripts", \
                           required=True)
    argparser.add_argument("-os", metavar="out_filtered_single", type=str, help="Output FASTA file with the filtered transcripts w/ single peak", \
                           required=True)
    argparser.add_argument("-pk", metavar="pickle_file", type=str, help="Output pickle file", required=True)    

    return argparser.parse_args()


def filter_transcripts(ref_genome, gff_file, peak_file):
    chroms = {}
    for seqr in SeqIO.parse(ref_genome, 'fasta'):
        if mito_gff_RE.search(seqr.description) is not None:
            chroms['chrMito'] = seqr.seq
            continue

        chr_name = chr_gff_RE.search(seqr.description).groups()[0]
        chr_name = 'chr%02d' % roman2int[chr_name]

        chroms[chr_name] = seqr.seq
        
    sgd_gff = HTSeq.GFF_Reader(gff_file)
    all_transcripts = []
    genes, nc_genes, genes_filtered = [], [], []
    features = defaultdict(lambda: defaultdict(list))

    feat_set = set() # TMP

    for feat in sgd_gff:
        if feat.attr.get('Parent') is None:
            feat_set.add(feat.type)
        # feat_set.add((feat.type, feat.attr.get('Parent')))
        if 'Parent' in feat.attr:
            features[feat.attr['Parent']][feat.type].append(feat)
        if feat.type == 'gene':
            if chr_RE.match(feat.iv.chrom):
                feat.iv.chrom = 'chr%02d' % roman2int[feat.iv.chrom[3:]]

            genes.append(feat)

        elif feat.type in [ "ncRNA", "rRNA", "snRNA", "snoRNA", "tRNA", "transposable_element_gene"]:
            if chr_RE.match(feat.iv.chrom):
                feat.iv.chrom = 'chr%02d' % roman2int[feat.iv.chrom[3:]]

            genes.append(feat)

    for gene in genes: 
        if gene.iv.chrom == '2-micron': # No sequence for 2-mu in the reference
            continue 

        if gene.type in [ 'gene', 'transposable_element_gene', 'pseudogene' ]:
            exons = [ chroms[gene.iv.chrom][e.iv.start:e.iv.end] \
                      for e in features[gene.name]['CDS'] ]
            ## Record the genes with no introns
            if not (len(features[gene.name]['intron']) or \
               len(features[gene.name]['five_prime_UTR_intron'])):

                genes_filtered.append(gene)
        else:
            exons = [ chroms[gene.iv.chrom][e.iv.start:e.iv.end] \
                      for e in features[gene.name]['noncoding_exon'] ]

        # for f in features[gene.name].keys():
        #     if f not in [ 'noncoding_exon', 'CDS' ]:
        #         print gene.name, gene.type, f
            
        ## Behold!
        seq = reduce(operator.add, exons, Seq.Seq(""))
            
        if gene.iv.strand == '-':
            seq = seq.reverse_complement()

        seqr = SeqRecord.SeqRecord(seq, id=gene.name, description="")
        all_transcripts.append(seqr)


    transcripts = [ g.name for g in genes_filtered ]

    all_peaks = {}
    peak_array = HTSeq.GenomicArrayOfSets("auto")
    peak_gff = HTSeq.GFF_Reader(peak_file)

    for peak in peak_gff:
        all_peaks[peak.name] = peak
        if (int(peak.attr['tfill_ypd_rep1']) >= 10 and int(peak.attr['tfill_ypd_rep2']) >= 10 and int(peak.attr['tfill_ypd_rep3']) >= 10):
            peak_array[peak.iv] = peak.name

    transcripts_single = []
    peak_info = {}
    for gene in genes_filtered:
        ## Filter out only peaks from tuples returned by steps()
        peaks = [ p for p in peak_array[gene.iv].steps() \
                  if (type(p[-1]) is str and p[-1].startswith('peak_'))]

        ## Record the peak information, including within-transcript coordinates, for pickling
        peak_info[gene.name] = []
        for peak in peaks:
            peak = all_peaks[peak[-1]]
            pi = ( peak.name, tuple(int(peak.attr['tfill_ypd_rep%i'%i]) for i in range(1, 3+1)), \
              (gene.iv.chrom, peak.iv.start - gene.iv.start, peak.iv.end - gene.iv.start, peak.iv.strand) )

            peak_info[gene.name].append(pi)

        n_pk = len(peaks)

        if n_pk == 1:
            transcripts_single.append(gene.name)

    return peak_info, all_transcripts, transcripts, transcripts_single
    

if __name__ == '__main__':
    args = parse_args()
    peak_info, transcripts_full, transcripts_filtered, transcripts_single_peak = filter_transcripts(args.r, args.t, args.p)

    pickle.dump(peak_info, file(args.pk, 'wb'))
    SeqIO.write(transcripts_full, file(args.ot, 'w'), 'fasta')
    file(args.o, 'w').writelines((g + '\n' for g in transcripts_filtered))
    file(args.os, 'w').writelines((g + '\n' for g in transcripts_single_peak))
    # SeqIO.write(transcripts_filtered, file(args.o, 'w'), 'fasta')
