#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
#


# This script generates a png file for each entry in a VCF file. The file
# displays the reads in a bam file around the VCF entry (+-WINDOW_SIZE).
# Reallocation of memry works now
from cmath import e
import os
import sys
import subprocess
import argparse
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
from matplotlib import gridspec

from vaf_from_pysam import vaf_from_pileup

# Optional dependencies
try:
	import pysam
except ModuleNotFoundError:
	# Handle the case when pysam is not installed
	print("pysam is not installed. Please install it for VAF annotation.")

try:
	import pysamstats
except ModuleNotFoundError:
	# Handle the case when pysamstats is not installed
	print("pysamstats is not installed. Please install it for VAF annotation.")


def get_args():
	argument_parser = argparse.ArgumentParser(
		description='This script generates a png file for each entry in a vcf file, a bed file or a manually specified region.' )
	argument_parser.add_argument('--bam_paths', metavar='STR', type=str,
			default=None, required=True, help='input list of bam files separated by comma. Maximum 3 BAM files')
	argument_parser.add_argument('--bam_names', metavar='STR', type=str,
			default=None, required=True, help='input list of names separated by comma. Same length as BAM files')
	argument_parser.add_argument('--ref', metavar='FILE', type=str,
			required=True, help='input reference genome file (fastq format)')
	argument_parser.add_argument('--exclude_flag', metavar='INT', type=str,
			default='3840', required=False, help='Exclude the reads with corresponding SAM flags, [default = %(default)s]')
	argument_parser.add_argument('--map_quality', metavar='INT', type=int,
			default='20', required=False, help='Minimum mapping quality for the reads, [default = %(default)s]')
	argument_parser.add_argument('--base_quality', metavar='INT', type=int,
			default='13', required=False, help='Minimum base quality for the variant, [default = %(default)s]')
	argument_parser.add_argument('--max_depth_plot', metavar='INT', type=int,
			default='500', required=False, help='Maximum read depth used to plot the high coverage region, [default = %(default)s]')
	argument_parser.add_argument('--vaf', action='store_true',
			help='Include the VAF of the central position in the plot title. Requires reference genome')
	argument_parser.add_argument('--for_gSmVs', action='store_true',
			help='For the gSmVs workflow used internally in DKFZ, the VAFs are directly sourced from the input VCF.')
	argument_parser.add_argument('--vcf', metavar='FILE', type=str,
			default=None, help='input vcf file ( as an alternative use --bed )')
	argument_parser.add_argument('--bed', metavar='FILE', type=str,
			default=None, help='input bed file ( as an alternative use --vcf )')
	argument_parser.add_argument('--annotations', metavar='FILE', type=str,
			default=None, help='Annotation track in bed format is indexed with a tabix. The fourth column could contain the annotation text for the segments. A comma can separate multiple files.')
	argument_parser.add_argument('--annotation_names', metavar='STR', type=str,
			default=None, help='annotation names separated by comma. Same length as annotation files')
	argument_parser.add_argument('--prefix', metavar='PREFIX', type=str,
			default="smvplot", help='target directory and file name prefix for generated output files, [default = %(default)s]')
	argument_parser.add_argument('--window', metavar='N', type=int,
			default=100, help='the output file for position X will show the region [X-window,X+window], [default = %(default)d]')
	argument_parser.add_argument('--samtoolsbin', metavar='N', type=str,
			default="samtools", help='the path to the samtools binary, [default = %(default)s]')
	argument_parser.add_argument('--tabixbin', metavar='N', type=str,
			default="tabix", help='the path to the tabix binary, [default = %(default)s]')
	argument_parser.add_argument('region', nargs='?', type=str,
			default=None, help='syntax either \'chr:start-end\' or \'chr:center\', use --vcf or --bed for more convenience')
	argument_parser.add_argument('--plot_dir', metavar='DIR', type=str,
			default="./", help='subfolder for the plots')
	argument_parser.add_argument('--out_format', metavar='STR', type=str,
			default='pdf', required=False, help='Output format of the plot, [default = %(default)s]')
	argument_parser.add_argument('--out_filename', metavar='STR', type=str,
			default=None, required=False, help='Output filename of the plot, [default = %(default)s]')
	argument_parser.add_argument('--ref_base', metavar='STR', type=str,
			default='', required=False, help='Reference base for the variant entry, [default = %(default)s]')
	argument_parser.add_argument('--alt_base', metavar='STR', type=str,
			default='', required=False, help='Alternate base for the variant entry, [default = %(default)s]')
	argument_parser.add_argument('--sort_by_variant', action='store_true',
			help='Sort the reads based on the variant base')


	parsed_arguments = argument_parser.parse_args()

	if len(parsed_arguments.bam_paths.split(',')) != len(parsed_arguments.bam_names.split(',')):
		print("-"*80)
		print("ERROR: Length of '--bam_paths' should match the length of '--bam_names'")
		print("-"*80)
		argument_parser.print_help()
		sys.exit()

	# check if the BAM files exit in the system, else remove the names from bam_names
	bam_paths = parsed_arguments.bam_paths.split(',')
	bam_names = parsed_arguments.bam_names.split(',')

	bam_paths_new, bam_names_new = [], []

	for path, name in zip(bam_paths, bam_names):
			if path[-3:] == 'bam' or path[-4:] == 'cram':
					bam_paths_new.append(path)
					bam_names_new.append(name)
			else:
					print(f"WARNING: BAM file {path} does not end with '.bam' or '.cram'")

	if len(bam_paths_new) == 0:
			print("-"*80)
			print("ERROR: No BAM files found")
			print("-"*80)
			argument_parser.print_help()
			sys.exit()
	parsed_arguments.bam_paths = ','.join(bam_paths_new)
	parsed_arguments.bam_names = ','.join(bam_names_new)

	if len(parsed_arguments.annotations.split(',')) != len(parsed_arguments.annotation_names.split(',')):
		print("-"*80)
		print("ERROR: Length of '--annotations' should match the length of '--annotation_names'")
		print("-"*80)
		argument_parser.print_help()
		sys.exit()

	return(parsed_arguments)

class ReferenceBuffer(object):

	def __init__( self, filename, chromosome ):

		self.samtools_call = [ parsed_arguments.samtoolsbin, "faidx", filename ]
		self.chromosome    = chromosome;
		self.offset        = 0
		self.sequence      = ""

	def __getitem__( self, pos ):

		if len(self.sequence) > pos - self.offset >= 0:

			return self.sequence[ pos - self.offset ]

		else:

			# To accomodate the smaller contigs
			try:
				region = "%s:%i-%i" % ( self.chromosome, pos - 1000, pos + 1000 )
				call   = self.samtools_call + [ region ]
				output = subprocess.check_output( call, encoding='UTF-8' )
				self.sequence = "".join( output.split('\n')[1:] ).upper()
				self.offset   = max( 0, pos - 1000 )
			except:
					region = self.chromosome
					call   = self.samtools_call + [ region ]
					output = subprocess.check_output( call, encoding='UTF-8' )
					self.sequence = "".join( output.split('\n')[1:] ).upper()
					self.offset = 1

			if pos - self.offset >= len(self.sequence):
				return "N"
			else:
				return self.sequence[ pos - self.offset ]


def variation_af(bam, chr, pos, map_quality, base_quality):
	"""Calculate the allele frequency of a variant in a set of BAM files.

	Args:
		bams (list): List of BAM files.
		chr (str): Chromosome.
		pos (int): Position.

	Returns:
		float: Allele frequency.
	"""

	# Iterate over the BAM files

	pos_zero = pos - 1

	variant_info = pysamstats.load_variation(bam,
					  chrom = chr,
						start = pos_zero,
						end = pos,
						fafile = parsed_arguments.ref,
						min_mapq = map_quality,
						min_baseq = base_quality)

	# variant info details from pysamstats out
	# REF to https://github.com/alimanfoo/pysamstats/blob/2e0980933494d9ce71639eed8c739ce9c9aa4617/pysamstats/opt.pyx#L515

	# Iterate over the list and tuple[1] == pos_zero
	try:
		variant_info_zero = [v for v in variant_info if v[1] == pos_zero][0]
	except:
		variant_info_zero = None
		Warning("No read support at position %s:%s in the %s BAM file" % (chr, pos, bam))

	read_count_ref = 0
	read_count_alt = 0

	if variant_info_zero is not None:
		# ref coverage count
		read_count_ref = variant_info_zero[3]

		if variant_info_zero[7] > 0 or variant_info_zero[11] > 0:
			# If it a mismatch
			if variant_info_zero[7] > variant_info_zero[11]:
				read_count_alt = variant_info_zero[7]
			# If it is an insertion
			elif variant_info_zero[11] > variant_info_zero[7]:
				read_count_alt = variant_info_zero[11]
		# A deletion perhaps
		else:
			# Check the next position
			try:
				variant_info_one = [v for v in variant_info if v[1] == pos][0]
			except:
				variant_info_one = None
				Warning("No read support at position %s:%s in the %s BAM file" % (chr, pos, bam))

			if variant_info_one is not None:
				if variant_info_one[9] > 0:
					read_count_ref = variant_info_one[3]
					read_count_alt = variant_info_one[9]

	variant_af = float(read_count_alt) / (read_count_ref) if (read_count_ref) > 0 else 0.0
	return variant_af


def get_annotations(region):
	annotations = []
	if parsed_arguments.annotations:
		annotation_files = parsed_arguments.annotations.split(',')
		for annotation_file in annotation_files:
			call = [parsed_arguments.tabixbin, annotation_file, region]
			output = subprocess.check_output(call, encoding='UTF-8')

			if not output:
				if region[:3] == "chr":
					call[2] = call[2][3:]
				else:
					call[2] = "chr" + call[2]

				output = subprocess.check_output(call, encoding='UTF-8')

			annotations.append(
				[[line.split('\t')[3], int(line.split('\t')[1]), int(line.split('\t')[2])]
				 for line in output.split('\n') if len(line.split('\t')) >= 3]
			)

	return annotations if annotations else None



def parse_cigar( cigar, pos, r_left, r_right):

	cigar_struct = [["", None]]
	# Adding a flanking region
	r_left = r_left - 10
	r_right = r_right + 10

	for char in cigar:
		# If the char is a number(in str form), add it to the last number as string
		# If the char is a letter, then convert the last number to int and add the letter
		# This cycle completes lenght of the cigar and the type of the cigar
		if "0" <= char <= "9":

			if type(cigar_struct[-1][0]) != str:

				cigar_struct.append([ "", None ])

			cigar_struct[-1][0] = cigar_struct[-1][0] + char

		else:

			cigar_struct[-1][0] = int( cigar_struct[-1][0] )
			cigar_struct[-1][1] = char

	if cigar_struct[-1][1] == None:

		del cigar_struct[-1]

	abs_pos = pos
	rel_pos = 0

	# If the first cigar is soft clipped, then the absolute position is shifted
	if cigar_struct[0][1] == 'S':
		abs_pos -= cigar_struct[0][0]

	cigar_struct2 = []

	# Iterate over the cigar structure
	for n,t in cigar_struct:

		# If the cigar is a match or a soft clipped, then add the position to the list
		if t in ['M','S']:

			for m in range(n):
				if r_left <= abs_pos+m <= r_right:
					cigar_struct2.append( (t,abs_pos+m, rel_pos+m) )

			abs_pos += n
			rel_pos += n

		elif t == 'D':

			for m in range(n):
				if r_left <= abs_pos+m <= r_right:
					cigar_struct2.append( (t,abs_pos+m, rel_pos) )

			abs_pos += n

		elif t == 'I':

			cigar_struct2.append( (t,abs_pos, rel_pos+m) )

			for m in range(1,n):
				if r_left <= abs_pos <= r_right:
					cigar_struct2.append( ('i',abs_pos, rel_pos+m) )

			rel_pos += n

		elif t == "N":
			for m in range(n):
				if r_left <= abs_pos+m <= r_right:
					cigar_struct2.append((t, abs_pos+m, rel_pos))
			abs_pos += n

	return cigar_struct2



def plot_histogram( cigars, ax, ax_id):

	if ax_id == 0:
		hist_color_codes = ['lightblue', 'darkblue', 'blue']
	elif ax_id == 1:
		hist_color_codes = ['lightcoral', 'darkred', 'red']
	elif ax_id == 2:
		hist_color_codes = ['lightgreen', 'darkgreen', 'green']

	points_deleted   = []
	points_clipped   = []
	points_unclipped = []

	for cigar in cigars:

		for (c_type, abs_pos, rel_pos) in cigar:

			if c_type in ['M']:
				points_deleted.append( abs_pos )
				points_clipped.append( abs_pos )
				points_unclipped.append( abs_pos )
			elif c_type in 'D':
				points_clipped.append( abs_pos )
				points_unclipped.append( abs_pos )
			elif c_type in ['S']:
				points_unclipped.append( abs_pos )

	points_deleted = sorted( points_deleted )
	points_clipped = sorted( points_clipped )
	points_unclipped = sorted( points_unclipped )

	min_unclipped = min( points_unclipped ) if len(points_unclipped) > 0 else 0
	max_unclipped = max( points_unclipped ) if len(points_unclipped) > 0 else 0

	histogram_deleted = [0] * ( 1 + max_unclipped - min_unclipped )
	histogram_clipped = [0] * ( 1 + max_unclipped - min_unclipped )
	histogram_unclipped = [0] * ( 1 + max_unclipped - min_unclipped )

	for p in points_deleted:

		histogram_deleted[ p - min_unclipped ] += 1

	for p in points_clipped:

		histogram_clipped[ p - min_unclipped ] += 1

	for p in points_unclipped:

		histogram_unclipped[ p - min_unclipped ] += 1


	original_x         = range( min_unclipped, max_unclipped + 1 )
	subdivided_x       = [ a+b for a in original_x for b in [-0.33,0.33]]

	smooth_unclipped = [ x for x in histogram_unclipped for y in [0,1] ]
	smooth_clipped   = [ x for x in histogram_clipped for y in [0,1] ]
	smooth_deleted   = [ x for x in histogram_deleted for y in [0,1] ]

	ax.fill_between( subdivided_x, 0, smooth_unclipped, color=hist_color_codes[0] )
	ax.fill_between( subdivided_x, 0, smooth_clipped, color=hist_color_codes[1] )
	ax.fill_between( subdivided_x, 0, smooth_deleted, color=hist_color_codes[2] )
	ax.set_ylim( ymin=0 )



def plot_cigars(cigars, sequences, reverses, ax, reference_function, region_center, region_chrom, read_names, read_names_vaf):
	patches = []

	right_limits = [0] * len(cigars)
	basepair_colors = { 'A':"#009600", 'C':"#3030fe", 'G':"#d17105", 'T':"#ff0000", 'N':"#00ffff" }
	rname_inside_region = []

	if parsed_arguments.sort_by_variant:
		# Sort cigars based on non-reference variants at the region center
		# The cigar list has the following structure: [ (cigar_type, absolute_position, relative_position), ... ]
		# The cigar_type can be 'M' (match), 'I' (insertion), 'D' (deletion), 'S' (soft clip), 'N' (skipped region)
		# The absolute_position is the position in the reference genome
		# The relative_position is the position in the read sequence
		# The cigar list is sorted based on the non-reference variants at the region center
		# The region center is the position of the variant in the reference genome
		# The reference_function is the reference genome sequence
		# Now separated reads based on the non-reference variants at the region center to group similar reads together
		# Create two lists: one for the reads that have the non-reference variant at the region center and one for the reads that do not have the non-reference variant at the region center
		
		cigar_with_variants = []
		cigar_without_variants = []
		seq_order_with_variants = []
		seq_order_without_variants = []
		rname_with_variants = []
		rname_without_variants = []

		for i, (read_cigar, sequence, rname) in enumerate(zip(cigars, sequences, read_names)):
			# Check if the region center is in the read
			for position_cigar in read_cigar:
				c_type, abs_pos, rel_pos = position_cigar
				if abs_pos == region_center:
					rname_inside_region.append(rname)
					if (c_type == 'M' and sequence[rel_pos] != reference_function[abs_pos]) and sequence[rel_pos] == parsed_arguments.alt_base:
						cigar_with_variants.append(read_cigar)
						seq_order_with_variants.append(i)
						rname_with_variants.append(rname)
						break
				elif abs_pos == (region_center + 1) and (len(parsed_arguments.alt_base) > 1 or len(parsed_arguments.ref_base) > 1):
					rname_inside_region.append(rname)
					if c_type == 'I' or c_type == 'D':
						cigar_with_variants.append(read_cigar)
						seq_order_with_variants.append(i)
						rname_with_variants.append(rname)
						break
			
			else:
				cigar_without_variants.append(read_cigar)
				seq_order_without_variants.append(i)
				rname_without_variants.append(rname)
		
		
		# Sort the reads with the non-reference variant at the region center based on the non-reference variant
		cigars_sorted = cigar_with_variants + cigar_without_variants
		sequences_sorted = [sequences[i] for i in seq_order_with_variants] + [sequences[i] for i in seq_order_without_variants]
		read_names_sored = rname_with_variants + rname_without_variants
	else:
		cigars_sorted = cigars
		sequences_sorted = sequences
		read_names_sored = read_names

	for cigar, sequence, reverse, rname in zip(cigars_sorted, sequences_sorted, reverses, read_names_sored):

		for j in range(len(cigars)):
			if right_limits[j] <= cigar[0][1]:
				line = j
				break

		clipped_cigar = [ (t,a,r) for (t,a,r) in cigar if t != 'S' ]

		right_limits[line] = cigar[-1][1] + 10

		# If the read is not in the region, then color it grey
		# If the read is in the region and used for VAF calculation, then color it grey
		# rest of the reads are colored white
		if rname in read_names_vaf or rname not in rname_inside_region:
			read_color = "#c8c8c8"
		else:
			read_color = "#f0f0f0"


		ax.barh( -line, clipped_cigar[-1][1]-clipped_cigar[0][1]+1, height=1, left=clipped_cigar[0][1]-0.5, color=read_color, linewidth=0 )

		for (c_type,abs_pos,rel_pos) in cigar:
			if c_type == 'D':
				ax.barh( -line, 1, height=1, left=abs_pos-0.5, color="#505050", linewidth=0 )
			elif c_type == 'M' and sequence[rel_pos] != reference_function[abs_pos]:
				ax.barh( -line, 1, height=1, left=abs_pos-0.5, color=basepair_colors[sequence[rel_pos]], linewidth=0, alpha=0.5 )
			elif c_type == 'S' and sequence[rel_pos] != reference_function[abs_pos]:
				ax.barh( -line, 1, height=1, left=abs_pos-0.5, color=basepair_colors[sequence[rel_pos]], linewidth=0, alpha=0.1 )
			elif c_type == 'N':
				ax.barh( -line, 1, height=1, left=abs_pos-0.5, color="#f0fdff", linewidth=0 )


		for (c_type,abs_pos,rel_pos) in cigar:
			if c_type == 'I':
				ax.barh( -line, 0.6, height=1, left=abs_pos-0.8, fill=False, linewidth=1 )

		ax.barh( -line, cigar[-1][1]-cigar[0][1]+1, height=1, left=cigar[0][1]-0.5, fill=False, linewidth=0 ) # linewidth changed to 0

		if reverse:
			patches.append( matplotlib.patches.Polygon( [(cigar[0][1]-0.5,-line),(cigar[0][1]-0.5,1-line),(cigar[0][1]-1.5,0.5-line)] ) )
			patches.append( matplotlib.patches.Polygon( [(cigar[-1][1]+0.5,1-line),(cigar[-1][1]+1,1-line),(cigar[-1][1]+0.5,0.5-line)] ) )
			patches.append( matplotlib.patches.Polygon( [(cigar[-1][1]+0.5,-line),(cigar[-1][1]+1,-line),(cigar[-1][1]+0.5,0.5-line)] ) )
		else:
			patches.append( matplotlib.patches.Polygon( [(cigar[-1][1]+0.5,-line),(cigar[-1][1]+0.5,1-line),(cigar[-1][1]+1.5,0.5-line)] ) )
			patches.append( matplotlib.patches.Polygon( [(cigar[0][1]-0.5,1-line),(cigar[0][1]-1,1-line),(cigar[0][1]-0.5,0.5-line)] ) )
			patches.append( matplotlib.patches.Polygon( [(cigar[0][1]-0.5,-line),(cigar[0][1]-1,-line),(cigar[0][1]-0.5,0.5-line)] ) )

	collection = matplotlib.collections.PatchCollection( patches, linewidths=0.5, edgecolors="black", facecolors="yellow" )
	ax.add_collection( collection )

	ax.set_ylim( ymin=( 1 - len([ l for l in right_limits if l > 0 ]) ), ymax=1 )
	ax.set_yticks([])



def plot_region(region_chrom, region_center, region_left, region_right, plot_title, VAFs=[]):
	region_string = "%s:%i-%i" % (region_chrom, region_left, region_right)
	annotations = get_annotations(region_string)
	basepair_colors = {'A': "#009600", 'C': "#3030fe", 'G': "#d17105", 'T': "#ff0000", 'N': "#00ffff"}

	bam_paths = parsed_arguments.bam_paths.split(",")
	bam_names = parsed_arguments.bam_names.split(",")

	annotation_idx = []

	if len(bam_paths) == 3:
		fig = plot.figure(figsize=(19.2, 16.2 + len(annotations)))
		rows = 10
		cols = 3
		h_ratios = [1, 5, 15, 2, 5, 15, 2, 5, 15]
		ax_indices = [1, 2, 4, 5, 7, 8]
		if annotations:
			ann_index = 10
			for i in range(0, len(annotations)):
				ax_indices += [ann_index + i * 2]
				annotation_idx.append(ann_index + i * 2)
	elif len(bam_paths) == 2:
		fig = plot.figure(figsize=(19.2, 10.8 + len(annotations)))
		rows = 7
		cols = 3
		h_ratios = [1, 5, 15, 2, 5, 15]
		ax_indices = [1, 2, 4, 5]
		if annotations:
			ann_index = 7
			for i in range(0, len(annotations)):
				ax_indices += [ann_index + i * 2]
				annotation_idx.append(ann_index + i * 2)
	else:
		fig = plot.figure(figsize=(19.2, 5.4 + len(annotations)))
		rows = 4
		cols = 3
		h_ratios = [1, 5, 15]
		ax_indices = [1, 2]
		if annotations:
			ann_index = 4
			for i in range(0, len(annotations)):
				ax_indices += [ann_index + i]
				annotation_idx.append(ann_index + i)

	if annotations:
		rows += 2 * len(annotations)
		h_ratios += [2, 1] * len(annotations)

	h_ratios += [1]

	grid = gridspec.GridSpec(rows, cols, height_ratios=h_ratios, hspace=0,
							 width_ratios=[1, 42, 1], wspace=0,
							 left=0, right=1, bottom=0, top=1)

	ax = [plot.subplot(grid[i, 1]) for i in ax_indices]

	reference_buffer = ReferenceBuffer(parsed_arguments.ref, region_chrom)
	visible_basepairs = [reference_buffer[i] for i in range(region_left, region_right + 1)]

	for idx, bam in enumerate(bam_paths):
		samtools_call = (parsed_arguments.samtoolsbin, "view",
						 "-q", str(parsed_arguments.map_quality),
						 "-F", parsed_arguments.exclude_flag, bam, region_string)

		samtools_output = subprocess.check_output(samtools_call, encoding='UTF-8')
		samtools_output = [line.split('\t') for line in samtools_output.split('\n')]
		samtools_reads = [line for line in samtools_output if len(line) > 5]

		region_read_counts = len(samtools_reads)
		if region_read_counts > parsed_arguments.max_depth_plot:
			print("*" * 80)
			print("Warning: too many reads in region %s, %s reads" % (region_string, region_read_counts))
			print("Warning: Taking a random of %s reads for plotting and histogram is still based on the original counts" % parsed_arguments.max_depth_plot)
			print("*" * 80)
			random_indices = random.sample(range(region_read_counts), parsed_arguments.max_depth_plot)
			samtools_reads_r = [samtools_reads[i] for i in sorted(random_indices)]

			samtools_read_plot = samtools_reads_r
			samtools_read_hist = samtools_reads
		else:
			samtools_read_plot = samtools_reads
			samtools_read_hist = samtools_reads


		plot_title_vaf = plot_title
		# Add all the read names to a list
		read_names_vaf = [ read[0] for read in samtools_read_plot ]

		if parsed_arguments.vaf:
			if parsed_arguments.for_gSmVs:
				if VAFs[idx] == ".":
					VAFs[idx] = "0.0"
				else:
					VAFs[idx] = str(round(float(VAFs[idx]), 4))
				plot_title_vaf = "%s - VAF=%s" % (plot_title, VAFs[idx])
			else:
				#vaf = variation_af(bam, region_chrom, region_center, parsed_arguments.map_quality, parsed_arguments.base_quality)
				# vaf_from_pileup(args.cram_file, args.ref_file, args.chrom, args.pos, args.var_ref, args.var_alt)
				read_names_vaf, alt_count, total_depth, vaf = vaf_from_pileup(bam,
																			parsed_arguments.ref,
																			region_chrom,
																			region_center,
																			parsed_arguments.ref_base,
																			parsed_arguments.alt_base,
																			parsed_arguments.map_quality,
																			parsed_arguments.base_quality)
				vaf = str(round(vaf, 4))
				plot_title_vaf = "%s - VAF=%s (%s/%s)" % (plot_title, vaf, alt_count, total_depth)

		# Plot the histogram and the cigars
		plot_histogram([parse_cigar(read[5], int(read[3]), region_left, region_right) for read in samtools_read_hist if read[5] != "*"], 
						ax[idx * 2],
						idx)
		plot_cigars([parse_cigar(read[5], int(read[3]), region_left, region_right) for read in samtools_read_plot if read[5] != "*"],
					[read[9] for read in samtools_read_plot if read[5] != "*"],
					[bool(int(read[1]) & 0x10) for read in samtools_read_plot if read[5] != "*"],
					ax[idx * 2 + 1], reference_buffer, region_center, region_chrom,
					[read[0] for read in samtools_read_plot if read[5] != "*"],
					read_names_vaf)


		if len(bam_paths) >= 2:
			if idx == 0:
				ax[idx * 2].set_title("%s - %s" % (plot_title_vaf, bam_names[0]))
			elif idx == 1:
				ax[idx * 2].set_title("%s - %s" % (plot_title_vaf, bam_names[1]))
			elif idx == 2:
				ax[idx * 2].set_title("%s - %s" % (plot_title_vaf, bam_names[2]))
		else:
			ax[idx * 2].set_title("%s - %s" % (plot_title_vaf, bam_names[0]))
		ax[idx * 2].set_xticks([])
		ax[idx * 2].yaxis.set_tick_params(labelleft=True, labelright=True)
		ax[idx * 2].ticklabel_format(style='plain', axis='x', useOffset=False)

		ax[idx * 2 + 1].xaxis.set_tick_params(width=0)
		ax[idx * 2 + 1].set_xticks([i for i in range(region_left, region_right + 1)])
		ax[idx * 2 + 1].xaxis.set_ticklabels(visible_basepairs, fontfamily="monospace")

		for tick in ax[idx * 2 + 1].get_xticklabels():
			tick.set_color(basepair_colors[tick._text])

		for axis in ax[idx * 2:(idx + 1) * 2]:
			axis.axvline(region_center - 0.5, color="black", linewidth=0.5)
			axis.axvline(region_center + 0.5, color="black", linewidth=0.5)
			for x in range(10, parsed_arguments.window, 10):
				axis.axvline(region_center - x - 0.5, color="black", linewidth=0.25)
				axis.axvline(region_center + x + 0.5, color="black", linewidth=0.25)

	for axis in ax:
		axis.set_xlim(xmin=region_left - 0.5, xmax=region_right + 0.5)

	if annotations:

		for i, anns in enumerate(annotations):

			ax_index = -(len(annotations)) + i

			ax_title = parsed_arguments.annotation_names.split(',')[i]

			ax[ax_index].set_title( "Annotations from " + ax_title, loc="left")
			ax[ax_index].set_xticks([])
			ax[ax_index].set_yticks([])
			ax[ax_index].axis('off')

			for ann in anns:
				ann[1] = max(ann[1], region_left)
				ann[2] = min(ann[2], region_right)

				ax[ax_index].barh(0, ann[2] - ann[1] + 1, height=1, left=ann[1] - 0.5, color="#c8c8c8", edgecolor="black", linewidth=0.5)
				ax[ax_index].text(float(ann[1] + ann[2]) / 2.0, 0.1, ann[0], ha='center', va='center')

			# Add a horizontal line to separate annotations
			#ax[ax_index].axhline(y = -1, color='black', linewidth=0.5)


def main():

	global parsed_arguments
	parsed_arguments = get_args()

	os.makedirs(parsed_arguments.plot_dir, exist_ok = True)

	if parsed_arguments.out_filename:
		out_filename = "%s.%s" %(parsed_arguments.out_filename, parsed_arguments.out_format)
	else:
		out_filename = "%s_%s_%i.%s" %(parsed_arguments.prefix, region_chrom, region_center, parsed_arguments.out_format)


	if parsed_arguments.region:

		region_chrom = ':'.join(parsed_arguments.region.split(':')[:-1])

		if len(parsed_arguments.region.split(':')[-1].split('-')) == 2:

			region_left   = int(parsed_arguments.region.split(':')[-1].split('-')[0])
			region_right  = int(parsed_arguments.region.split(':')[-1].split('-')[1])
			region_center = ( region_left + region_right ) // 2

		else:

			region_center = int(parsed_arguments.region.split(':')[-1])
			region_left   = region_center - parsed_arguments.window
			region_right  = region_center + parsed_arguments.window

		plot_title = "%s:%s" % ( region_chrom, region_center )
		if parsed_arguments.ref_base and parsed_arguments.alt_base:
			plot_title += " %s>%s" % (parsed_arguments.ref_base, parsed_arguments.alt_base)

		# For variants with multiple SNVs, split them and assume the first base 
		# in the reference and alternate alleles are the different and 
		# take is as the ref and alt base for sorting
		if len(parsed_arguments.ref_base) == len(parsed_arguments.alt_base):
			parsed_arguments.ref_base = parsed_arguments.ref_base[0]
			parsed_arguments.alt_base = parsed_arguments.alt_base[0]

		#print(region_chrom, region_left, region_center, region_right)

		plot_region( region_chrom, region_center, region_left, region_right, plot_title )

		plot_output_file = os.path.join(parsed_arguments.plot_dir, out_filename)
		plot.savefig( plot_output_file )
		plot.clf()
		plot.cla()
		plot.close()



	if parsed_arguments.vcf:

		vcf_columns = {}

		for line in open(parsed_arguments.vcf, 'r' ):
			line = line.rstrip('\n')

			if line[:1] == "#":
				for i,col in enumerate( line[1:].rstrip('\n').split('\t') ):
					vcf_columns[col] = i

			elif line[:1] != "#" and "CHROM" in vcf_columns and "POS" in vcf_columns:
				region_chrom  = line.split('\t')[ vcf_columns["CHROM"] ]
				region_center = int( line.split('\t')[ vcf_columns["POS"] ] )
				region_left   = region_center - parsed_arguments.window
				region_right  = region_center + parsed_arguments.window
				region_ref    = line.split('\t')[ vcf_columns["REF"] ]
				region_alt    = line.split('\t')[ vcf_columns["ALT"] ]

				plot_title = "%s:%s" % ( region_chrom, region_center )
				plot_title += " %s>%s" % (region_ref, region_alt)

				# For variants with multiple SNVs, split them and assume the first base 
				# in the reference and alternate alleles are the different and 
				# take is as the ref and alt base for sorting
				if len(region_ref) == len(region_alt):
					parsed_arguments.ref_base = region_ref[0]
					parsed_arguments.alt_base = region_alt[0]

				parsed_arguments.ref_base = region_ref
				parsed_arguments.alt_base = region_alt

				if "VEP_Most_Severe_Consequence" in vcf_columns and line.split('\t')[ vcf_columns["VEP_Most_Severe_Consequence"] ] != '.':
					plot_title += " %s" % line.split('\t')[ vcf_columns["VEP_Most_Severe_Consequence"] ].replace("_", " ").rstrip('\n')            # .rstrip('\n') added to avoid new line in title
				if "HUGO_Symbol" in vcf_columns and line.split('\t')[ vcf_columns["HUGO_Symbol"] ] != '.':
					plot_title += " in %s" % line.split('\t')[ vcf_columns["HUGO_Symbol"] ].rstrip('\n')		# .rstrip('\n') added to avoid new line in title

				VAFs = []
				if parsed_arguments.for_gSmVs:
					VAFs = [
						line.split('\t')[ vcf_columns["Control_VAF"] ],
						line.split('\t')[ vcf_columns["Tumor_VAF"] ],
						line.split('\t')[ vcf_columns["RNA_VARIANT_AF"] ],
					]

				plot_region( region_chrom, region_center, region_left, region_right, plot_title, VAFs)
				plot_output_file = os.path.join(parsed_arguments.plot_dir, out_filename)
				plot.savefig( plot_output_file )
				plot.clf()
				plot.cla()
				plot.close()



	if parsed_arguments.bed:

		for line in open(parsed_arguments.bed, 'r' ):

			if line[:1] != "#":

				region_chrom  = line.split('\t')[0]
				region_left   = int(line.split('\t')[1])
				region_right  = int(line.split('\t')[2])
				region_center = ( region_left + region_right ) // 2

				plot_title = "%s:%s" % ( region_chrom, region_center )

				plot_region( region_chrom, region_center, region_left, region_right, plot_title )
				plot_output_file = os.path.join(parsed_arguments.plot_dir, out_filename)
				plot.savefig( plot_output_file )
				plot.clf()
				plot.cla()
				plot.close()

if __name__ == "__main__":
	main()