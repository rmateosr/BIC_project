#!/usr/bin/env python3
# ABOUTME: Bypass modkit for nanopolish-era samples — convert the per-read methylation TSV
# ABOUTME: + phased BAM into the stage-5 merged-by-chr TSV that linebylineextractingreadinfo_HG002.R reads.
#
# Output path (matches coordfixer.py output):
#   <output_folder>/modkit_referenced_splitbam_mergedbychr/modkit_extract_<chr>_merged_filtered_coordmod.tsv
#
# Output columns (tab-separated, no header — matches coordfixer.py):
#   read_id   cpg_position   chrom   methylation_status(0/1)   haplotype(H1|H2|noH)
#
# Rows are sorted by read_id (grouping so R's single-pass buffer works), then by position
# within read. Ambiguous calls (code 'x') are dropped, matching nanopolish's confidence threshold.
#
# The raw nanopolish TSV (tabix-indexed, 7 columns) is assumed to be an aggregated per-read
# summary format with column 5 encoding CpG calls as "<offset><code>..." where offset is bp
# from the previous CpG start (first offset is 0) and code is m=methylated, u=unmethylated,
# x=ambiguous. CpG genomic position = row.start + cumulative_offset.
import sys, os, re, subprocess, shutil

if len(sys.argv) != 5:
    sys.stderr.write(
        "Usage: nanopolish_to_mergedbychr.py <chr> <meth_tsv_gz> <phased_bam> <output_folder>\n"
    )
    sys.exit(2)

CHROM, METH_TSV, BAM, OUT = sys.argv[1:]
OUT_DIR = os.path.join(OUT, "modkit_referenced_splitbam_mergedbychr")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PATH = os.path.join(OUT_DIR, f"modkit_extract_{CHROM}_merged_filtered_coordmod.tsv")

# ------------------------------------------------------------
# 1. Build read_id -> haplotype map from the phased BAM (primary + supplementary kept out;
#    -F 260 skips unmapped + secondary; HP tag is consistent across a read's alignments anyway).
# ------------------------------------------------------------
sys.stderr.write(f"[{CHROM}] scanning BAM for HP tags...\n")
hap = {}
HP_RE = re.compile(r"^HP:i:([12])$")
with subprocess.Popen(
    ["samtools", "view", "-F", "260", BAM, CHROM],
    stdout=subprocess.PIPE,
    text=True,
) as p:
    assert p.stdout is not None
    for line in p.stdout:
        f = line.split("\t")
        qname = f[0]
        if qname in hap:
            continue
        h = "noH"
        for t in f[11:]:
            m = HP_RE.match(t.rstrip())
            if m:
                h = "H" + m.group(1)
                break
        hap[qname] = h
    rc = p.wait()
if rc != 0:
    sys.stderr.write(f"[{CHROM}] samtools view exit={rc}\n")
    sys.exit(rc)
sys.stderr.write(f"[{CHROM}] BAM reads mapped to chr: {len(hap)}\n")

# ------------------------------------------------------------
# 2. Stream tabix-extracted rows and emit one row per CpG to an unsorted temp.
# ------------------------------------------------------------
tmp_unsorted = OUT_PATH + ".unsorted.tmp"
tmp_sorted = OUT_PATH + ".sorted.tmp"
OFFSET_RE = re.compile(r"(\d+)([muxMUX])")

n_in = n_out = n_x = n_noread = 0
with open(tmp_unsorted, "w") as fout, subprocess.Popen(
    ["tabix", METH_TSV, CHROM], stdout=subprocess.PIPE, text=True
) as p:
    assert p.stdout is not None
    for line in p.stdout:
        n_in += 1
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        c, start_s, _end_s, rid, calls_s = parts[:5]
        if c != CHROM:
            continue
        try:
            start = int(start_s)
        except ValueError:
            continue
        h = hap.get(rid)
        if h is None:
            n_noread += 1
            h = "noH"  # TSV had a read not in BAM — treat as unphased, keep the data
        pos = start
        for m in OFFSET_RE.finditer(calls_s):
            d = int(m.group(1))
            code = m.group(2).lower()
            pos += d  # first match has d=0 so pos stays at start
            if code == "m":
                prob = "1"
            elif code == "u":
                prob = "0"
            else:
                n_x += 1
                continue
            fout.write(f"{rid}\t{pos}\t{CHROM}\t{prob}\t{h}\n")
            n_out += 1
    rc = p.wait()
if rc != 0:
    sys.stderr.write(f"[{CHROM}] tabix exit={rc}\n")
    sys.exit(rc)
sys.stderr.write(
    f"[{CHROM}] tabix rows: {n_in}; cpg emitted: {n_out}; ambiguous dropped: {n_x}; "
    f"reads not in BAM (kept as noH): {n_noread}\n"
)

# ------------------------------------------------------------
# 3. Sort by read_id (grouping so R's single-pass buffer works), then by position.
#    Use GNU sort — handles >memory input via on-disk merge.
# ------------------------------------------------------------
sys.stderr.write(f"[{CHROM}] sorting...\n")
with open(tmp_sorted, "w") as fout:
    subprocess.check_call(
        ["sort", "-t", "\t", "-k1,1", "-k2,2n", "-S", "4G", tmp_unsorted],
        stdout=fout,
    )
os.remove(tmp_unsorted)

# ------------------------------------------------------------
# 4. Atomic rename to final path. No header (matches coordfixer.py behavior — R's
#    linebylineextractingreadinfo_HG002.R silently discards line 1 of its input).
# ------------------------------------------------------------
os.rename(tmp_sorted, OUT_PATH)
sys.stderr.write(f"[{CHROM}] wrote {OUT_PATH}\n")
