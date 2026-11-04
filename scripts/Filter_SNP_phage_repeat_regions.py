#!/usr/bin/env python3
import argparse
from multiprocessing import Pool, cpu_count
from typing import List, Tuple, Optional, Set

def parse_args() -> argparse.Namespace:
	"""
	Parse command‑line arguments.
	"""
	parser = argparse.ArgumentParser(
		description='Filter VCF to exclude SNPs in repeat and phage regions, using multiprocessing.'
	)
	parser.add_argument('--vcf', type=str, required=True, help='Input VCF file to be filtered.')
	parser.add_argument('--position', type=str, required=True, help='TSV of start/end positions (tab‑delimited) of phage & repeat regions.')
	parser.add_argument('--output', type=str, default='Filtered_NoPhage_NoRepeat.vcf', help='Path to write filtered VCF (default: %(default)s).')
	parser.add_argument('--threads', type=int, default=8, help='Number of processes to use for filtering (default: 8).')
	return parser.parse_args()

def load_regions(position_file: str) -> List[Tuple[int, int]]:
	"""
	Read a two‑column TSV of start/end positions and return a list of (start, end) tuples.
	"""
	regions: List[Tuple[int,int]] = []
	with open(position_file) as f:
		for line in f:
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			start_str, end_str = line.split('\t')[:2]
			regions.append((int(start_str), int(end_str)))
	return regions

def record_in_regions(pos: int, regions: List[Tuple[int, int]]) -> bool:
	"""
	Return True if pos falls within any (start, end) region.
	"""
	return any(start <= pos <= end for start, end in regions)

def _filter_line(args: Tuple[str, List[Tuple[int,int]]]) -> Optional[str]:
	"""
	Helper for Pool.map: returns the line if it's NOT in any region, else None.
	"""
	line, regions = args
	pos = int(line.split('\t')[1])
	if not record_in_regions(pos, regions):
		return line
	return None

def filter_vcf_parallel(
	vcf_file: str,
	regions: List[Tuple[int, int]],
	processes: int
) -> Tuple[List[str], List[str]]:
	"""
	Read VCF, collect headers and data lines, then parallel‑filter data lines.
	Returns (headers, included_data).
	"""
	headers: List[str] = []
	data_lines: List[str] = []
	with open(vcf_file) as f:
		for raw in f:
			line = raw.rstrip('\n')
			if line.startswith('#'):
				headers.append(line)
			else:
				data_lines.append(line)

	# build the arg‑tuples for mapping
	args_iter = ((ln, regions) for ln in data_lines)
	with Pool(processes=processes) as pool:
		results = pool.map(_filter_line, args_iter)

	# drop None entries
	included = [ln for ln in results if ln is not None]
	return headers, included

def write_output(
	output_file: str,
	headers: List[str],
	data_lines: List[str]
) -> None:
	"""
	Write headers and data_lines to output_file.
	"""
	with open(output_file, 'w') as out:
		for h in headers:
			out.write(h + '\n')
		for ln in data_lines:
			out.write(ln + '\n')

def main():
	args = parse_args()
	regions = load_regions(args.position)

	if args.threads > 1:
		headers, included = filter_vcf_parallel(
			args.vcf, regions, args.threads
		)
	else:
		# If only 1 thread requested, fall back to serial version:
		headers, included = filter_vcf_parallel(
			args.vcf, regions, 1
		)

	write_output(args.output, headers, included)
	print(f"Total SNP Count: {len(included)}")

if __name__ == '__main__':
	main()
