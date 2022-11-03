#!/usr/bin/env python3

# Irenaeus Chan
# 04/08/2022
# Looks for evidence of complex reads within BAM of a sample
# e.g. Given two independent variants: chr2 25240411 AG>A, chr2 25240418 G>T within the same
# Sample, we want to know if there is evidence that both these variants come from the SAME
# read/clone. This would indicate a complex variant that wasn't caught by the caller.

# Specifications
# Both variants have to appear on a single read
# If there is evidence of both variants appearing in two different reads, they are two separate
# clones and shouldn't be considered a complex variant

# Steps:
# 1. Figure out what the combined variant should "look" like within the BAM file
# e.g. chr2 25240411 AG>A, chr2 25240418 G>T
# 11 12 13 14 15 16 17 18
#  A  G  -  -  -  -  -  G
#  A  -  -  -  -  -  T          <- Due to the deletion, the T is moved over

# More Complicated
# chr17 60663124 GT>G
# chr17 60663125 TA>G
# chr17 60663009 TC>T
# chr17 60663009 TCCATGGCCAAG>T
# chr17 60663136 TCAAAAGATCCAGAACCACTTGAAG>T
# chr17 60663143 AT>A
# chr17 60663144 TC>T
# chr17 60663182 CT>C
# chr17 60663182 CTT>C

# 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85
#  G  T  A  -  -  -  -  -  -  -  -  -  T  C  A  A  A  A  G  A  T  C  C  A  G  A  A  C  C  A  C  T  T  G  A  A  G  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  C  T  T

# Several Possibilities
# 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85
#  G  G  -  -  -  -  -  -  -  -  -  T  C  A  A  A  A  G  A  T  C  C  A  G  A  A  C  C  A  C  T  T  G  A  A  G  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  C  T  T

import pysam
import os
import sys
import argparse
import string
import re

 #SNPs are easy... it's just same position, single Nucleotide change
 #INDELs are hard because it involves shifting left (Deletion) or right (Insertion)...
 #This means we should probably classify what the mutation is, so we are aware of what
 #we need to do with it. Aka, if it's Insertion, we will have to shift the sequence right (+1)

class Variant(object):
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = str(chrom)
        self.pos = int(pos)
        self.ref = str(ref)
        self.alt = str(alt)
        if self.is_snp(): self.mut_type = "snp"
        if self.is_insertion(): self.mut_type = "insertion"
        if self.is_deletion(): self.mut_type = "deletion"
        if self.is_complex(): self.mut_type = "complex"

    @classmethod
    def createVariant(self, variant):
        return Variant(variant.split(":")[0], variant.split(":")[1], variant.split(":")[2], variant.split(":")[3])

    def is_snp(self):
        if len(self.ref) == 1 and len(self.alt) == 1: return True
        else: return False

    def is_indel(self):
        if len(self.ref) | len(self.alt) > 1: return True
        else: return False

    def is_insertion(self):
        if (len(self.ref) == 1 and len(self.alt) > 1 and self.ref[0] == self.alt[0]): return True
        else: return False

    def is_deletion(self):
        if (len(self.ref) > 1 and len(self.alt) == 1 and self.ref[0] == self.alt[0]): return True
        else: return False

    def is_complex(self):
        if (len(self.ref) > 1 and len(self.alt) > 1): return True
        elif (self.is_indel() and self.ref[0] != self.alt[0]): return True
        else: return False

    def toString(self):
        return(self.chrom + ":" + str(self.pos) + ":" + self.ref + ":" + self.alt)

def cigarToCigarettes(cigarstring): return(re.findall(r'(\d+)([A-Z]{1})', cigarstring))

def doesRecordContainVariant(record, variant):
    if record.is_unmapped == True: return False                                         # If not mapped, can have other issues
    if record.reference_name != variant.chrom: return False                             # Ensure that we only look at the current chromosome of interest
    if record.reference_start > variant.pos-1: return False                             # Is our variant within the read itself (-1 because reference_start is based on index 0)
    if record.reference_end < variant.pos + len(variant.alt): return False              # Is our variant within the read (+ ALT because if this is an INDEL we need the INDEL to be in the read to check)
    return True

# Iterates through the CIGAR string and determines if the variant position needs to be altered due to INDEL prior
# to the variant itself within the record. If an insertion occurs before the variant, we would have to increase by
# the size of the insertion. Similarly, if a deletion occurs, we would have to decrease by the size of the deletion.
# e.g.
#18M2I55M54H chr17 60663227 60663247 20 AAC:A
#TGTGCCTACTAATTCAAC  AAACACTGTCATGGACCAAAAAAATNTGAAGATGTCAACTCCTGGCCAAATGAAA    REF
#TGTGCCTACTAATTCAACATAAACACTGTCATGGACCAAAAAAATNTGAAGATGTCAACTCCTGGCCAAATGAAA    RECORD
#                     ^
#32M1D46M2D26M26H chr17 60663168 60663247 79 AAC:A
#GCTAAAGCCCTGACTTTAAGGATACATGATTCTTTGAATAATAGCCTTCCAATTGGCCTTGTGCCTACTAATTCAACAAACACTGTCATGGACCAAAAAAATTTGAA  REF
#GCTAAAGCCCTGACTTTAAGGATACATGATTC TTGAATAATAGCCTTCCAATTGGCCTTGTGCCTACTAATTCAACAA  ACTGTCATGGACCAAAAAAATTTGAA  RECORD
def calculatePosInRecord(reference_pos, reference_start, cigarstring):
    variant_pos = reference_pos - reference_start
    leftover_bases = variant_pos
    for bases, cigar_code in cigarToCigarettes(cigarstring):
        if cigar_code in ['H', 'S']: continue
        if cigar_code == "I":
            variant_pos = variant_pos + int(bases)
        if cigar_code == "D":
            variant_pos = variant_pos - int(bases)
        if leftover_bases - int(bases) <= 0: return variant_pos
        else:
            leftover_bases = leftover_bases - int(bases)
            continue

#15M4D114M chr17 60663052 60663068 12 TCTCA:C TTCGTAGCAATGCCTC

# Iterates through the CIGAR string to determine whether the INDEL variant belongs to the INDELs found within the record.
# If a record has two or three deletions, we need to make sure one of the deletions is the specific deletion we are interested in.
# Insertion is fairly simple because we just make sure whether or not the denoted insertion matches our insertion of interest.
# Deletion is harder since we need to know what replaced the nucleotides we lost. So we need to check the reference to generate
# a replacement string which we know the record should be at the specific variant position.
# e.g. Deletion
# AACACTGTCAT
# A  ACTGTCAT     The AC is replaced with another AC. So we need to make sure at the location of the deletion, the record has AAC and not AACAC
def checkIndelWithCigar(indel_pos_in_record, variant, record, reference):
        leftover_bases = indel_pos_in_record
        indel_start = -1
        found_variant = False
        for bases, cigar_code in cigarToCigarettes(record.cigarstring):
            if cigar_code in ['H', 'S']: continue
            if cigar_code == "I" and variant.is_insertion() and len(variant.alt)-1 == int(bases) and indel_pos_in_record == indel_start+1:
                #print(record.cigarstring + ' ' + str(indel_pos_in_record) + ' ' + str(indel_start) + ' ' + record.query_alignment_sequence + ' ' + record.query_alignment_sequence[indel_start:indel_start+len(variant.alt)])
                #if (record.query_alignment_sequence[indel_start:indel_start+len(variant.alt)] == variant.alt):
                found_variant = True
            if cigar_code == "D" and variant.is_deletion() and len(variant.ref)-1 == int(bases) and indel_pos_in_record == indel_start+1:
                #ref_one = pysam.faidx(reference, variant.chrom+":"+str(variant.pos)+"-"+str(variant.pos))
                #ref_two = pysam.faidx(reference, variant.chrom+":"+str(variant.pos+(len(variant.ref)))+"-"+str(variant.pos+2*(len(variant.ref)-1)))
                #deletion = ref_one.split("\n")[1]+ref_two.split("\n")[1]
                #print(record.cigarstring + ' ' + str(indel_pos_in_record) + ' ' + str(indel_start) + ' ' + record.query_alignment_sequence + ' ' + record.query_alignment_sequence[indel_start:indel_start+len(variant.ref)])
                #if(record.query_alignment_sequence[indel_start:indel_start+len(variant.ref)] == deletion):
                found_variant = True
            if leftover_bases - int(bases) < 0: return found_variant
            else:
                leftover_bases = leftover_bases - int(bases)
                indel_start = indel_start + int(bases)
                continue

def findVariant(bam, variant, reference):
    record_list = []
    if variant.is_complex():
        print(pysam.faidx(reference, variant.chrom+":"+str(variant.pos-1)+"-"+str(variant.pos)).split("\n")[1])
        #if len(variant.alt) == 1
    for record in bam.fetch(variant.chrom, variant.pos-500, variant.pos+500):
        if doesRecordContainVariant(record, variant):
            variant_start_pos = calculatePosInRecord(variant.pos, record.reference_start, record.cigarstring)
            if variant.is_snp() and record.query_alignment_sequence[variant_start_pos-1] == variant.alt: record_list.append(record)
            if variant.is_complex():
                print(record.cigarstring + record.reference_name + str(record.reference_start) + " " + str(variant.pos) + " " + str(variant_start_pos) + " " + variant.ref + ":" + variant.alt + " " + record.query_alignment_sequence)
                #chr17:76737015:CGG:TGA    8036706596_S82
                #chr17:60663068:TCTCA:C    8036708641_S127
            if variant.is_insertion() and "I" in record.cigarstring:
                if checkIndelWithCigar(variant_start_pos, variant, record, reference):
                    record_list.append(record)
            if variant.is_deletion() and "D" in record.cigarstring:
                if checkIndelWithCigar(variant_start_pos, variant, record, reference):
                    record_list.append(record)
    return(record_list)

def variantToBAM(bam, variant, reference, outfile, complex):
    with pysam.AlignmentFile(outfile, "wb", header=bam.header) as outf:
        if complex != None:
            set_before = set()
            for variant in complex:
                variant = Variant.createVariant(variant)
                list_of_records = findVariant(bamfile, variant, reference)
                if set_before == set():
                    set_before = set(list_of_records)
                set_before = set_before & set(list_of_records)
            for record in list(set_before):
                outf.write(record)
        else:
            for record in findVariant(bam, Variant.createVariant(variant), reference):
                outf.write(record)

def variantCount(bam, variant, reference, complex):
    if complex != None:
        print(complex)
        dict_of_complex_variants = {}
        set_before = set()
        for variant in complex:
            variant = Variant.createVariant(variant)
            if variant.toString() not in dict_of_complex_variants:
                dict_of_complex_variants[variant.toString()] = list()
            list_of_records = findVariant(bamfile, variant, reference)
            dict_of_complex_variants[variant.toString()] = list_of_records
            if set_before == set():
                set_before = set(list_of_records)
            set_before = set_before & set(list_of_records)
        for variants, list_of_records in dict_of_complex_variants.items():
            print("Number of reads detected with", variants, ":",len(list_of_records))
        print("Number of reads detected with ALL variants:",len(list(set_before)))
    else:
        print("Number of reads detected with", variant, ":", len(findVariant(bam, Variant.createVariant(variant), reference)))

def readsToBAM(bam, variant, outfile, complex):
    with pysam.AlignmentFile(outfile, "wb", header=bam.header) as outf:
        if complex != None:
            print(complex)
            set_before = set()
            list_of_records = []
            for variant in complex:
                variant = Variant.createVariant(variant)
                for record in bam.fetch(variant.chrom, variant.pos-500, variant.pos+500):
                    if doesRecordContainVariant(record, variant):
                        list_of_records.append(record)
                if set_before == set():
                    set_before = set(list_of_records)
                set_before = set_before & set(list_of_records)
            for record in list(set_before):
                outf.write(record)
        else:
            variant = Variant.createVariant(variant)
            for record in bam.fetch(variant.chrom, variant.pos-500, variant.pos+500):
                if doesRecordContainVariant(record, variant):
                    print(record)
                    outf.write(record)


def parser_setup():
    parser = argparse.ArgumentParser(description='Pulls the reads that match specific variants of interest. Used as a replacement for checking in the IGV tool.\n \
                                                For example: chr2:25240418:G:T would pull all the reads in the BAM that have this specific variant.')
    parser.add_argument('function', metavar='TOOL', type=str, nargs='+',help='VariantToBAM, VariantCount')
    parser.add_argument('-r', '--reference', dest='reference', type =argparse.FileType('r'), help = 'Reference File', required = True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-v', '--variant', dest='variant', type = str, help = 'Variant in format: CHROM:POS:REF:ALT')
    group.add_argument('-c', '--complex', dest='complex', type = str, nargs='+', help = 'A list of variants that form a complex variant e.g. chr2:25240411:AG:A, chr2:25240418:G:T')
    parser.add_argument('-b', '--bam', dest='BAM', type=argparse.FileType('r'), help = 'The BAM/SAM/CRAM file to search for the variant', required = True)
    parser.add_argument('-o', '--outfile', dest='outfile', type = argparse.FileType('w'), help = 'The file to write results')

    return parser

if __name__ == "__main__":
    parser = parser_setup()
    args = parser.parse_args()

    ext = os.path.splitext(args.BAM.name)[-1].lower()
    if (ext == '.bam'):
        bamfile = pysam.AlignmentFile(args.BAM, "rb")
    elif (ext == '.sam'):
        bamfile = pysam.AlignmentFile(args.BAM, "r")
    elif (ext == '.cram'):
        bamfile = pysam.AlignmentFile(args.BAM, "rc")
    else:
        parser.print_help()
        exit(1)

    for tool in args.function:
        if tool == 'VariantToBAM':
            if args.outfile == None:
                print("ERROR: Output File Required\n")
                parser.print_help()
                exit(1)
            variantToBAM(bamfile, args.variant, args.reference.name, args.outfile.name, args.complex)
        elif tool == 'VariantCount':
            variantCount(bamfile, args.variant, args.reference.name, args.complex)
        elif tool == 'ReadsToBAM':
            if args.outfile == None:
                print("ERROR: Output File Required\n")
                parser.print_help()
                exit(1)
            readsToBAM(bamfile, args.variant, args.outfile.name, args.complex)


#chr4	105235980	NA	A	T
#chr17	60656798	NA	C	CTT
#chr17	60663247	NA	AAC	A

#chr17:60663068:TCTCA:C

#chr17:60663067:TTCTC:T
#chr17:60663072:A:C

#chr2:25240411:AG:A
#chr2:25240418:G:T

#chr2:25247149:CAA:C
#chr2:25247150:AAC:T

#chr17:60663077:AT:A
#chr17:60663078:TT:A

#139M4D12M chr1760662928 60663068 136 TCTCA:C CATGCATAGATTTGTTGAGTTCTGGGATAAATTTTTTCTTATTTGTTTTACCTTCTTATTTTTCAGTCACTGGAGGAGGATCCATGGCCAAGGGTGAATNCTAAGGACCATATACCTGCCCTGGTTNGTAGCAATGCCTCGAGAATTTTTT
#15M4D114M chr17 60663052 60663068 12 TCTCA:C TTCGTAGCAATGCCTCGAGAATTTTTTAGAGGTTTCAGCTGAGATAGCTCGAGAGAATGTCCAAGGTGTAGTCATACCCTCAAAAGATCCAGAACCACTTGAAGAAAATTGNGCTAAAGCCCTGACTTT
