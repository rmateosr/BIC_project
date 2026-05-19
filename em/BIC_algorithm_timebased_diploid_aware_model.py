# ABOUTME: Driver for time-based ASM detection under the diploid per-allele marginal constraint.
# ABOUTME: Identical to BIC_algorithm_timebased.py except for the EM module it imports.
import argparse
import glob
import os
import warnings
from multiprocessing import Pool

import EMfunctions_timebased_diploid_aware_model as em


def main():
    warnings.filterwarnings("ignore", category=FutureWarning, module="numpy")

    parser = argparse.ArgumentParser(
        description="Time-based ASM detection using EM + BIC on multi-time-period data."
    )
    parser.add_argument(
        'input_paths', nargs='+',
        help='Input directories, one per time period (in chronological order)'
    )
    parser.add_argument('--output', '-o', type=str, required=True, help='Output TSV path')
    parser.add_argument('--windowsize', type=int, default=10, help='CpG sites per window')
    parser.add_argument('--workers', type=int, default=10,
                        help='Parallel worker processes (default 10 matches BIC ASM production).')
    parser.add_argument('--glob-pattern', type=str, default='methylationfraction_*_.tsv',
                        help='Glob pattern for per-region TSVs produced by dividemethylationintosmallerregions_updated_HG002.R.')
    args = parser.parse_args()

    T = len(args.input_paths)

    # Find region files for each time period
    files_by_t = {}
    for t, dir_path in enumerate(args.input_paths, 1):
        files_by_t[t] = sorted(glob.glob(os.path.join(dir_path, args.glob_pattern)))
        print(f"Time period {t}: {len(files_by_t[t])} region files in {dir_path}")

    # Verify all time periods have the same number of region files
    n_regions = len(files_by_t[1])
    for t in range(2, T + 1):
        if len(files_by_t[t]) != n_regions:
            raise ValueError(
                f"Time period {t} has {len(files_by_t[t])} files, "
                f"but time period 1 has {n_regions}. Region counts must match."
            )

    # Build tasks: one dict of paths per region
    tasks = []
    for region_idx in range(n_regions):
        paths_for_region = {t: files_by_t[t][region_idx] for t in range(1, T + 1)}
        tasks.append((paths_for_region, args.windowsize))

    # Incremental append as each region finishes: if the job is killed mid-run,
    # *.tmp retains every already-completed region. Atomic rename on clean exit.
    tmp_output = args.output + ".tmp"
    header_written = False
    with Pool(processes=args.workers) as pool, open(tmp_output, 'w') as f:
        for df in pool.imap_unordered(_process_region, tasks):
            df.to_csv(f, sep='\t', index=False, header=not header_written)
            header_written = True
            f.flush()
    os.replace(tmp_output, args.output)
    print(f"Output written to {args.output}")


def _process_region(task):
    """Wrapper for multiprocessing (single-arg for imap_unordered)."""
    paths_by_t, windowsize = task
    return em.EMBIC_bin_path_timebased(paths_by_t, windowsize=windowsize)


if __name__ == "__main__":
    main()
