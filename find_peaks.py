#!python3

# Copyright © 2022, Nan Huang

# This is a Python port to the original find_peaks program
# at https://github.com/owenjm/find_peaks

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
# USA

import argparse
import datetime
import math
import os
import random
import sys
from io import TextIOWrapper
from typing import Optional

version = "1.0.1"

parser = argparse.ArgumentParser(
    description="Simple FDR random permutation peak caller",
    allow_abbrev=True,
)
parser.add_argument(
    "--n",
    type=int,
    default=100,
    help="Number of iterations",
)
parser.add_argument(
    "--fdr",
    type=float,
    default=0.01,
    help="False discovery rate value",
)
parser.add_argument(
    "--frac",
    type=int,
    default=0,
    help="Number of random fragments to consider per iteration",
)
parser.add_argument(
    "--min_count",
    type=int,
    default=2,
    help="Minimum number of fragments to consider as a peak",
)
parser.add_argument(
    "--min_quant",
    type=float,
    default=0.95,
    help="Minimum quantile for considering peaks",
)
parser.add_argument(
    "--step",
    type=float,
    default=0.01,
    help="Stepping for quantiles",
)
parser.add_argument(
    "--unified_peaks",
    type=str,
    default="max",
    choices=["max", "min"],
    help="Method for calling peak overlaps (two options):\n\t'min': call minimum overlapping peak area\n\t'max': call maximum overlap as peak",
)
parser.add_argument(
    "--seed",
    type=int,
    default=0,
    help="Random seed",
)
parser.add_argument(
    "files",
    type=str,
    nargs="+",
    help="Input files in bedgraph or GFF format",
)
args = parser.parse_args()
random.seed(args.seed)


def load_gff(fn) -> list[tuple[str, int, int, float]]:
    global args
    sys.stderr.write(f"Reading input file: {fn} ...\n")
    with open(fn, "r") as f:
        lines = f.readlines()
    sys.stderr.write(f"Read {len(lines)} lines\n")
    parsed_result: list[tuple[str, int, int, float]] = list()
    for line in lines:
        ll = line.strip().split("\t")
        # skip empty lines
        if len(ll) < 4:
            continue
        if len(ll) == 4:
            # bedgraph
            chrom, start, end, score = ll
        else:
            # GFF
            chrom = ll[0]
            start, end, score = ll[3:6]
        parsed_result.append(
            (chrom, int(start), int(end), float(score) if score != "NA" else 0)
        )
    sys.stderr.write("Sorting ...\n")
    parsed_result = sorted(parsed_result, key=lambda x: (x[0], x[1]))
    return parsed_result


def find_quant(probes: list[tuple[str, int, int, float]]):
    global args
    total_coverage = 0
    frags: list[float] = list()
    for (chrom, start, end, score) in probes:
        if not score:
            continue
        total_coverage += end - start
        frags.append(score)
    frags = sorted(frags)
    sys.stderr.write(f"Total coverage was {total_coverage} bp\n")
    quants = [
        (q * args.step, int(q * args.step * len(frags)) - 1)
        for q in range(math.ceil(args.min_quant / args.step), math.ceil(1 / args.step))
    ]
    for cut_off, score_idx in quants:
        sys.stdout.write(f"\tQuantile {cut_off:0.2f}: {frags[score_idx]:0.2f}\n")
    peakmins = [frags[score_idx] for (_, score_idx) in quants]
    return peakmins


def call_peaks_unified_redux(
    iter_num: int,
    probes: list[tuple[str, int, int, float]],
    peakmins: list[float],
    peaks: Optional[
        dict[float, list[tuple[str, int, int, float, float, int, int]]]
    ] = None,
    peak_count: Optional[dict[float, dict[int, int]]] = None,
    peak_count_real: Optional[dict[float, dict[int, int]]] = None,
    real: bool = False,
):
    global args

    if real:
        sys.stderr.write("Calling real peaks ...\r")
    else:
        sys.stderr.write(f"Iteration {iter_num+1}: [processing ...]\r")
    if not peak_count:
        peak_count = dict()
    if not peak_count_real:
        peak_count_real = dict()
    if not peaks:
        peaks = dict()
    pstart = pend = pscore = inpeak = count = 0
    tmp_peak = list()
    old_chrom = ""
    for pm in peakmins:
        peaks[pm] = peaks.get(pm, list())
        peak_count[pm] = peak_count.get(pm, dict())
        peak_count_real[pm] = peak_count_real.get(pm, dict())
        for chrom, start, end, score in filter(lambda x: x[3], probes):
            if real:
                if chrom != old_chrom:
                    # Next chromosome
                    # (Peaks can't carry over chromosomes, but we don't use this shortcut when randomly shuffling)
                    pstart = pend = pscore = inpeak = count = 0
                    tmp_peak = list()
            old_chrom = chrom if real else old_chrom
            if not inpeak:
                if score >= pm:
                    # record new peak
                    pstart = start
                    pend = end
                    pscore = score * (end - start) / 1000
                    inpeak = 1
                    count += 1
                    if real:
                        tmp_peak.append(score)
                else:
                    continue
            else:
                if score >= pm:
                    # still in peak
                    count += 1
                    # Fragment score to deal with scoring peaks made from uneven sized fragments
                    fragment_score = score * (end - start) / 1000
                    if real:
                        tmp_peak.append(score)
                    pscore += fragment_score
                    pend = end
                else:
                    # Out of a peak
                    if count >= args.min_count:
                        # record peak
                        if real:
                            peak_count_real[pm][count] = (
                                peak_count_real[pm].get(count, 0) + 1
                            )
                            mean_pscore = round(pscore / (pend - pstart) * 1000, 2)
                            peaks[pm].append(
                                (
                                    chrom,
                                    pstart,
                                    pend,
                                    mean_pscore,
                                    pscore,
                                    count,
                                    pend - pstart,
                                )
                            )
                        else:
                            peak_count[pm][count] = peak_count[pm].get(count, 0) + 1
                    # reset
                    pstart = pend = pscore = inpeak = count = 0
                    if real:
                        tmp_peak = list()
    if real:
        return peaks, peak_count_real
    return peaks, peak_count


def find_randomised_peaks(
    probes: list[tuple[str, int, int, float]], peakmins: list[float]
):
    global args
    peak_count = None
    sys.stdout.write("Duplicating ...\n")
    pbs = probes.copy()
    sys.stdout.write("Calling peaks on input file ...\n")
    for iter_num in range(args.n):
        sys.stdout.write(f"Iteration {iter_num+1}: [shuffling]\r")
        if args.frac:
            pbs = pbs[: args.frac]
        # The built-in shuffle uses the same algorithm (Fisher-Yates)
        # as the original Perl program
        random.shuffle(pbs)
        _, peak_count = call_peaks_unified_redux(
            iter_num, pbs, peakmins, None, peak_count
        )
    return peak_count


def calculate_regressions(
    probes: list[tuple[str, int, int, float]],
    peakmins: list[float],
    peak_count: dict[float, dict[int, int]],
    file_handle: TextIOWrapper,
):
    global args
    log_scores = dict()
    regression = dict()
    num = len(probes)
    for pm in peakmins:
        file_handle.write(f"Peak min = {pm}\n")
        for c in peak_count[pm].keys():
            peak_count_avg = peak_count[pm][c] / args.n if peak_count[pm][c] else 0
            if not peak_count_avg:
                continue
            if args.frac:
                peak_count_avg *= num / args.frac
            log_scores[pm] = log_scores.get(pm, dict())
            log_scores[pm][c] = math.log10(peak_count_avg)
            file_handle.write(f"Peak size: {c}\tCount: {peak_count_avg}\n")
        # calculate exponential decay rates
        # y = ax+b for log(y)
        sumx = sumy = sumxy = sumx2 = 0
        n = 0
        for c in peak_count[pm].keys():
            if not peak_count[pm][c]:
                continue
            n += 1
            sumx += c
            sumy += log_scores[pm][c]
            sumxy += c * log_scores[pm][c]
            sumx2 += c**2
        if n < 2:
            continue
        mean_x = sumx / n
        mean_y = sumy / n
        mean_x2 = sumx2 / n
        mean_xy = sumxy / n
        a = (mean_xy - mean_x * mean_y) / (mean_x2 - mean_x**2)
        b = mean_y - a * mean_x

        # store values
        regression[pm] = (a, b)
        file_handle.write(f"regression: log(y) = {a}(x) + {b}\n")

        for c in peak_count[pm].keys():
            if not peak_count[pm][c]:
                continue
            a, b = regression[pm]
            logval = a * c + b
            val = 10**logval
            file_handle.write(
                f"lin regress: {c}\t{log_scores[pm][c]}\t{logval}\t{val}\n"
            )
        file_handle.write("\n")
    return log_scores, regression


def calculate_fdr(
    peakmins: list[float],
    regression: dict[float, tuple[float, float]],
    peak_count_real: dict[float, dict[int, int]],
    file_handle: TextIOWrapper,
):
    global args
    fdr: dict[float, dict[int, float]] = dict()
    peak_fdr_cutoff: dict[float, int] = dict()
    for pm in peakmins:
        fdr[pm] = fdr.get(pm, dict())
        # get regression variables
        if not regression.get(pm):
            continue
        a, b = regression[pm]
        for c in peak_count_real[pm].keys():
            if not peak_count_real[pm][c]:
                continue
            expect = 10 ** (a * c + b)
            real_count = peak_count_real[pm][c]
            fdr_conservative = expect / real_count
            fdr[pm][c] = fdr_conservative
    # print FDR rates
    file_handle.write("\n")
    for pm in peakmins:
        file_handle.write(f"Peak min = {pm}\n")
        for c in sorted(fdr[pm].keys()):
            file_handle.write(
                f"Peak size: {c}\tCount: {peak_count_real[pm][c]}\tFDR: {fdr[pm][c]}\n"
            )
            if fdr[pm][c] < args.fdr and not peak_fdr_cutoff.get(pm, 0):
                peak_fdr_cutoff[pm] = c
        # clumsy hack to prevent errors
        peak_fdr_cutoff[pm] = peak_fdr_cutoff.get(pm, int(1e10))
        file_handle.write("\n")

    for pm in peakmins:
        file_handle.write(
            f"Peak min {pm}: peak cutoff size for alpha = {args.fdr} was {peak_fdr_cutoff[pm]}\n\n"
        )
    return fdr, peak_fdr_cutoff


def find_significant_peaks(
    peaks: dict[float, list[tuple[str, int, int, float, float, int, int]]],
    peakmins: list[float],
    fdr: dict[float, dict[int, float]],
    peak_fdr_cutoff: dict[float, int],
):
    global args
    sig_peaks: list[tuple[str, int, int, float, float, int, int, float]] = list()
    sys.stderr.write("Selecting significant peaks ...\n")
    for pm in peakmins:
        for chrom, pstart, pend, mean_pscore, pscore, count, size in peaks[pm]:
            if count >= peak_fdr_cutoff[pm]:
                sig_peaks.append(
                    (
                        chrom,
                        pstart,
                        pend,
                        mean_pscore,
                        pscore,
                        count,
                        size,
                        fdr[pm][count],
                    )
                )
    return sig_peaks


def make_unified_peaks(
    sig_peaks: list[tuple[str, int, int, float, float, int, int, float]],
    out_file: str,
):
    global args
    # Unify overlapping peaks, and make significant peaks file
    unified_peaks: list[tuple[str, str, str, int, int, float, str, str, str]] = list()
    total = len(sig_peaks)
    i = 0
    sys.stderr.write("Combining significant peaks ...\n")

    # unroll chromosomes for speed
    for chrom in set([x[0] for x in sig_peaks]):
        sig_peaks_in_chrom = filter(lambda x: x[0] == chrom, sig_peaks)
        unified_peaks_chr: list[
            tuple[str, str, str, int, int, float, str, str, str]
        ] = list()
        for (
            chrom,
            start,
            end,
            score,
            total_score,
            count,
            peaklen,
            fdr,
        ) in sig_peaks_in_chrom:
            i += 1
            if i % 100 == 0:
                sys.stderr.write(f"{i/total*100:0.2f}% processed ...\r")
            # next if unified_peaks_chr already overlaps
            if any(x[3] < end and start < x[4] for x in unified_peaks_chr):
                continue
            # Grab all elements that overlap
            overlap = filter(lambda x: x[1] < end and start < x[2], sig_peaks_in_chrom)

            for (
                chrom_o,
                start_o,
                end_o,
                score_o,
                total_score_o,
                count_o,
                peaklen_o,
                fdr_o,
            ) in overlap:
                if args.unified_peaks == "min":
                    start = max(start, start_o)
                    end = min(end, end_o)
                else:
                    start = min(start, start_o)
                    end = max(end, end_o)
                score = max(score, score_o)
                fdr = min(fdr, fdr_o)

            unified_peaks_chr.append(
                (chrom, ".", ".", start, end, score, ".", ".", f"FDR={fdr}")
            )
        unified_peaks += unified_peaks_chr
    sys.stderr.write("Sorting unified peaks ...\n")
    unified_peaks = sorted(unified_peaks, key=lambda x: (x[0], x[3]))
    sys.stderr.write("Writing unified peaks file ...\n")
    with open(out_file, "w") as file_handle:
        for peak in unified_peaks:
            file_handle.write("\t".join([str(x) for x in peak]) + "\n")
    return unified_peaks


def main():
    global args
    date = datetime.datetime.now().strftime("%Y-%m-%d.%H-%M-%S")
    for fn in args.files:
        # path/to/file.bedgraph => path/to, file, .bedgraph
        dir, (name, ext) = (os.path.dirname(fn), os.path.splitext(os.path.basename(fn)))

        # Output file names
        fn_base_date = "peak_analysis." + name + "." + date
        base_dir = os.path.join(dir, fn_base_date)
        out = os.path.join(base_dir, name + f"-FDR{args.fdr:.2f}")
        out_peak_unified_track = out + ".peaks.gff"
        out_peaks = out + "_FDR-data"

        # Load gff/bedgraph data files
        probes = load_gff(fn)
        # Calculate min values of peaks for given quantiles
        peakmins = find_quant(probes)
        peak_count = find_randomised_peaks(probes, peakmins)
        assert peak_count
        # Make directory
        os.makedirs(base_dir, exist_ok=True)
        # Open peaks file for writing
        with open(out_peaks, "w") as f:
            # Write header
            f.write(f"FDR peak call v{version}\n\n")
            f.write(f"Input file: {fn}\n")
            log_scores, regression = calculate_regressions(
                probes, peakmins, peak_count, f
            )
            # peaks were only recorded if they were real
            peaks, peak_count_real = call_peaks_unified_redux(
                1, probes, peakmins, real=True
            )
            fdr, peak_fdr_cutoff = calculate_fdr(
                peakmins, regression, peak_count_real, f
            )
            sig_peaks = find_significant_peaks(peaks, peakmins, fdr, peak_fdr_cutoff)
            unified_peaks = make_unified_peaks(sig_peaks, out_peak_unified_track)
            sys.stderr.write(f"{len(unified_peaks)} peaks found.\n\n")


if __name__ == "__main__":
    main()
