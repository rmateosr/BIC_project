# ABOUTME: Builds per-chromosome union of CpG coordinates across two samples.
# ABOUTME: Output is one coord per line — consumed by dividemethylationintosmallerregions_updated_HG002.R.
"""
build_union_coords.py

Reads two read_format/${chr}_reads.tsv files (one per time point) and writes the
union of CpG coordinates found in the `startcoord` column (comma-separated ints)
to OUTPUT_FILE. The R region-splitter consumes this as `allcoords` so that
region boundaries match across samples — mandatory for time-based ASM pairing.
"""
import argparse
import os
import sys


def load_coords(path):
    coords = set()
    with open(path, 'r') as fh:
        header = fh.readline()
        # Expect columns: readlabel, haplotype, chrom, startcoord, status
        cols = header.rstrip('\n').split('\t')
        try:
            idx = cols.index('startcoord')
        except ValueError:
            sys.exit(f"startcoord column missing in {path}: {cols}")
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= idx:
                continue
            for c in fields[idx].split(','):
                c = c.strip()
                if not c:
                    continue
                try:
                    coords.add(int(c))
                except ValueError:
                    pass
    return coords


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reads-file', action='append', required=True,
                    help='Per-sample read_format/${chr}_reads.tsv; repeat for each sample.')
    ap.add_argument('--output', required=True, help='Destination union-coords file.')
    args = ap.parse_args()

    union = set()
    for p in args.reads_file:
        if not os.path.isfile(p):
            sys.exit(f"Missing reads file: {p}")
        union |= load_coords(p)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    tmp = args.output + '.tmp'
    with open(tmp, 'w') as fh:
        for c in sorted(union):
            fh.write(f"{c}\n")
    os.replace(tmp, args.output)
    print(f"Wrote {len(union)} coords -> {args.output}")


if __name__ == '__main__':
    main()
