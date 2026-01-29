#!/usr/bin/env python3
# Calculate VAF from CRAM file using pysam

import argparse
import pysam

import pysam

def vaf_from_pileup(bam_file, ref_file, chrom, pos, var_ref, var_alt, min_mqual, min_bqual, ignore_overlaps=True, compute_baq=True, redo_baq=False):
    # Open the BAM file
    samfile = pysam.AlignmentFile(bam_file, "rb", reference_filename=ref_file)
    fastafile = pysam.FastaFile(ref_file)
    
    # Get the pileup at the variant position
    pileup = samfile.pileup(chrom, pos - 1, pos + len(var_ref),
                            stepper="samtools",
                            min_base_quality=min_bqual,
                            min_mapping_quality=min_mqual,
                            ignore_overlaps=ignore_overlaps,
                            compute_baq=compute_baq,
                            redo_baq=redo_baq,
                            fastafile=fastafile)
    
    if not pileup:
        print(f"No reads found at {chrom}:{pos}")
        return None
    
    ref_count = 0
    alt_count = 0
    total_depth = 0
    read_names = []
    
    for pileupcol in pileup:
        if pileupcol.pos + 1 == pos:
            for pileupread in pileupcol.pileups:

                # Skip if query position is None
                if pileupread.query_position is None:
                    continue

                read_names.append(pileupread.alignment.query_name)
                ref_base = pileupread.alignment.query_sequence[pileupread.query_position]
                base_quality = pileupread.alignment.query_qualities[pileupread.query_position]

                if len(var_ref) == 1 and len(var_alt) == 1:
                    if ref_base == var_ref:
                        ref_count += 1
                    elif ref_base == var_alt:
                        alt_count += 1
                elif len(var_ref) == len(var_alt):
                    # MNV handling
                    mnv_len = len(var_alt)
                    qp = pileupread.query_position
                    if qp + mnv_len <= len(pileupread.alignment.query_sequence):
                        read_substr = pileupread.alignment.query_sequence[qp:qp+mnv_len] 
                        if read_substr == var_alt:
                            alt_count += 1
                        elif read_substr == var_ref:
                            ref_count += 1
                else:
                    # Insertion or deletion
                    indel_indicator = pileupread.indel
                    ins_or_del = 'Ins' if len(var_alt) > len(var_ref) else 'Del'
                    alt_len = len(var_alt) - 1
                    ref_len = len(var_ref) - 1

                    if indel_indicator == 0:
                        ref_count += 1
                    elif ins_or_del == 'Ins' and indel_indicator == alt_len :
                        alt_count += 1
                    elif ins_or_del == 'Del' and indel_indicator == -ref_len:
                        alt_count += 1
                    else:
                        pass
                    #print(f"Unexpected base {ref_base} at {chrom}:{pos}")

                # Print the base quality scores here
                #print(f"Base quality: {base_quality}\t{ref_base}\t{var_ref}\t{var_alt}\t{ref_count}\t{alt_count}")
                total_depth += 1
    
    if total_depth == 0:
        print(f"No reads covering the variant at {chrom}:{pos}")
        return None, None, None, None
    
    vaf = alt_count / total_depth

    return read_names, alt_count, total_depth, vaf

def main():
    parser = argparse.ArgumentParser(description="Calculate VAF from CRAM file using pysam")
    parser.add_argument("--cram_file", type=str, help="Path to the CRAM file", dest = 'cram_file')
    parser.add_argument("--ref_file", type=str, help="Path to the reference file", dest = 'ref_file')
    parser.add_argument("--chrom", type=str, help="Chromosome name", dest = 'chrom')
    parser.add_argument("--pos", type=int, help="Position", dest = 'pos')
    parser.add_argument("--var_ref", type=str, help="Reference allele", dest = 'var_ref')
    parser.add_argument("--var_alt", type=str, help="Alternate allele", dest = 'var_alt')
    parser.add_argument("--mqual", type=int, help="Minimum mapping quality", default=30)
    parser.add_argument("--bqual", type=int, help="Minimum base quality", default=13)
    
    args = parser.parse_args()
    
    read_names, alt_count, total_depth, vaf = vaf_from_pileup(args.cram_file, 
                                        args.ref_file, 
                                        args.chrom, 
                                        args.pos, 
                                        args.var_ref,
                                        args.var_alt,
                                        args.mqual,
                                        args.bqual)
    
    if vaf is not None:
        print(f"VAF at {args.chrom}:{args.pos} is {vaf:.4f}")

if __name__ == "__main__":
    main()