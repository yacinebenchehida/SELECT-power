#!/usr/bin/env bash
# =====================================================================================================
# Runs the SLiM side of the SLiM comparison (validation/slim.md).
# Requires SLiM 5 (https://messerlab.org/slim/; tested with SLiM 5.2). Run from the repository root:
#     bash validation/slim_runs.sh            (about 20 minutes; the largest runs need about 3 GB of memory each)
# then compare with SPOT:
#     node validation/slim_compare.js         (writes validation/slim_results.tsv, slim_power.tsv, slim_fpr.tsv)
#     python3 validation/slim_figure.py       (writes validation/slim_check.svg)
#
# Each line of SCENARIOS below is one SLiM run: R = 2,000 replicate populations of one scenario.
# Columns: name  system  NF  NM  P0  S  H  [seed]
#   name   : label used for the output file validation/slim_out/<name>.tsv (and read by slim_compare.js)
#   system : diploid, haploid or haplodiploid
#   NF, NM : breeding females and males per generation (haploid: NF individuals; NM unused)
#   P0     : starting frequency of the selected allele
#   S, H   : selection coefficient and dominance (S = 0 gives the neutral runs used to set the threshold)
#   seed   : random seed (default 12345); the extra neutral runs (_b, _c, _d) only add replicates for a more
#            precise check of the false-positive rate (see slim.md)
# All runs share the design of the default page: sampling at generations 0, 20, 40 and 90 (one generation
# per year, 1935/1955/1975/2025), 20/20/20/40 females per sampling date, perfect genotypes.
# =====================================================================================================
set -euo pipefail
SLIM=${SLIM:-slim}              # path to the SLiM executable (set SLIM=/path/to/slim if it is not on the PATH)
OUT=validation/slim_out
R=2000
mkdir -p "$OUT"

SCENARIOS="
diploid_N1000_s0.02         diploid       500  500  0.1   0.02  0.5
diploid_N1000_neutral       diploid       500  500  0.1   0     0.5
diploid_N1000_s0.05_rec     diploid       500  500  0.1   0.05  0
diploid_N200_s0.05          diploid       100  100  0.1   0.05  0.5
diploid_N200_neutral        diploid       100  100  0.1   0     0.5
haploid_N1000_s0.02         haploid       1000 0    0.05  0.02  0.5
haploid_N1000_neutral       haploid       1000 0    0.05  0     0.5
haplodiploid_N1000_s0.02    haplodiploid  500  500  0.1   0.02  0.5
haplodiploid_N1000_neutral  haplodiploid  500  500  0.1   0     0.5
haplodiploid_N2000_s0.05    haplodiploid  1000 1000 0.01  0.05  0.5
haplodiploid_N2000_neutral  haplodiploid  1000 1000 0.01  0     0.5
diploid_N1000_neutral_b     diploid       500  500  0.1   0     0.5  2
diploid_N1000_neutral_c     diploid       500  500  0.1   0     0.5  3
diploid_N1000_neutral_d     diploid       500  500  0.1   0     0.5  4
diploid_N200_neutral_b      diploid       100  100  0.1   0     0.5  2
diploid_N200_neutral_c      diploid       100  100  0.1   0     0.5  3
diploid_N200_neutral_d      diploid       100  100  0.1   0     0.5  4
haplodiploid_N1000_neutral_b haplodiploid 500  500  0.1   0     0.5  2
haplodiploid_N1000_neutral_c haplodiploid 500  500  0.1   0     0.5  3
haplodiploid_N1000_neutral_d haplodiploid 500  500  0.1   0     0.5  4
"

run_one() {   # arguments: name system NF NM P0 S H [seed]
	rm -f "$OUT/$1.tsv"
	"$SLIM" -s ${8:-12345} -d "SYSTEM='$2'" -d NF=$3 -d NM=$4 -d P0=$5 -d S=$6 -d H=$7 -d R=$R \
		-d "GENS=c(0,20,40,90)" -d "NSAMP=c(20,20,20,40)" -d "OUT='$OUT/$1.tsv'" \
		validation/slim_check.slim > /dev/null
	echo "done $1"
}
export -f run_one; export SLIM OUT R
# one run at a time by default (JOBS=2 bash validation/slim_runs.sh runs two at once if memory allows)
echo "$SCENARIOS" | grep -v '^$' | xargs -P ${JOBS:-1} -L 1 bash -c 'run_one "$@"' _
