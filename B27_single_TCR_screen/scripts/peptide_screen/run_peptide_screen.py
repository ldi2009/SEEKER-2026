#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Driver (peptide_screen): run the full Illumina peptide-only single-TCR screen
from config, sample by sample:

    1. extract peptide counts from the R2 clean FASTQ (01)
    2. merge R1 + R2 fullresult tables (02)
    3. annotate peptide origin against human/virus libraries (03)
    4. add pre-selection (R0) counts (04)
    5. keep only resolved-origin peptides (05)
    6. sort by R2 count and add synthesis primers (06)

Every command and its output is appended to a run log; a Markdown processing
report is written at the end. Progress and estimated remaining time are shown
after each step.

Usage:
    python run_peptide_screen.py --config configs/my_screen.yml
"""

import os
import sys
import gzip
import time
import shutil
import argparse
import subprocess
from datetime import datetime, timedelta

import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


class ProgressTracker:
    """Step-wise progress with average-based remaining-time estimates."""

    def __init__(self, total_steps, log_file):
        self.total_steps = total_steps
        self.current_step = 0
        self.start_time = time.time()
        self.log_file = log_file

    def start_step(self, step_name):
        self.current_step += 1
        self.step_start = time.time()
        msg = (f"\n{'='*60}\n"
               f"Step {self.current_step}/{self.total_steps}: {step_name}\n"
               f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
               f"{'='*60}")
        print(msg)
        self._log(msg)

    def end_step(self, step_name, details=""):
        elapsed = time.time() - self.step_start
        total_elapsed = time.time() - self.start_time
        avg = total_elapsed / self.current_step
        remaining = avg * (self.total_steps - self.current_step)
        msg = (f"\nStep {self.current_step}/{self.total_steps} done: {step_name}\n"
               f"This step: {self._fmt(elapsed)}\n"
               f"Total elapsed: {self._fmt(total_elapsed)}\n"
               f"Estimated remaining: {self._fmt(remaining)}\n")
        if details:
            msg += f"Details: {details}\n"
        msg += f"{'='*60}"
        print(msg)
        self._log(msg)

    @staticmethod
    def _fmt(seconds):
        if seconds < 60:
            return f"{seconds:.1f}s"
        if seconds < 3600:
            return f"{seconds/60:.1f}min"
        return f"{seconds/3600:.2f}h"

    def _log(self, msg):
        try:
            with open(self.log_file, 'a') as f:
                f.write(msg + '\n')
        except Exception:
            pass


def run_step(cmd, log_file):
    """Run a pipeline step, streaming its output into the log."""
    printable = ' '.join(str(c) for c in cmd)
    print(f"[CMD] {printable}")
    with open(log_file, 'a') as lf:
        lf.write(f"\n[CMD] {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        lf.write(f"{printable}\n")
        lf.write("-" * 40 + "\n")

    proc = subprocess.Popen(
        [str(c) for c in cmd],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    for line in proc.stdout:
        print(line, end='')
        with open(log_file, 'a') as lf:
            lf.write(line)
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"Step failed with exit code {proc.returncode}: {printable}")


def decompress_gz(gz_file, output_file):
    if os.path.exists(output_file):
        print(f"FASTQ already exists, skipping gunzip: {output_file}")
        return
    size_mb = os.path.getsize(gz_file) / (1024 * 1024)
    print(f"Decompressing {gz_file} ({size_mb:.1f} MB) -> {output_file}")
    with open(output_file, 'wb') as f_out:
        with gzip.open(gz_file, 'rb') as f_in:
            shutil.copyfileobj(f_in, f_out, length=1024 * 1024 * 8)


def run_sample(sample_name, sample_cfg, cfg, log_file, tracker):
    results_dir = cfg.get('output_dir', 'results')
    ref = cfg['reference']
    r0 = cfg.get('round0', {})
    max_reads = (cfg.get('extract', {}) or {}).get('max_reads')

    outputs = {}

    # ---- step 1: extract R2 ----
    tracker.start_step(f"[{sample_name}] extract peptide counts (R2)")
    r2_fullresult = sample_cfg.get('R2', {}).get('fullresult')
    if r2_fullresult and os.path.exists(r2_fullresult):
        print(f"R2 fullresult already provided, skipping extraction: {r2_fullresult}")
        outputs['r2_fullresult'] = r2_fullresult
    else:
        fq = sample_cfg['R2']['fastq']
        if fq.endswith('.gz'):
            fq_plain = os.path.splitext(fq)[0]
            decompress_gz(fq, fq_plain)
        else:
            fq_plain = fq
        out_name = f"{sample_name}-R2_plasmid_library.csv"
        cmd = [sys.executable, os.path.join(SCRIPT_DIR, '01_extract_peptide_counts.py'),
               fq_plain, out_name, results_dir]
        if max_reads:
            cmd += ['--max_reads', str(max_reads)]
        run_step(cmd, log_file)
        outputs['r2_fullresult'] = os.path.join(results_dir, f"fullresult_{out_name}")
    tracker.end_step("extract peptide counts (R2)")

    # ---- step 2: combine rounds ----
    tracker.start_step(f"[{sample_name}] combine R1 + R2")
    combined = os.path.join(results_dir, f"{sample_name}_R1R2_combined.csv")
    run_step([sys.executable, os.path.join(SCRIPT_DIR, '02_combine_rounds.py'),
              sample_cfg['R1']['fullresult'], outputs['r2_fullresult'],
              combined, sample_name], log_file)
    outputs['combined'] = combined
    tracker.end_step("combine R1 + R2")

    # ---- step 3: annotate origin ----
    tracker.start_step(f"[{sample_name}] annotate origin")
    annotated = os.path.join(results_dir, f"{sample_name}_R1R2_annotated.csv")
    cmd = [sys.executable, os.path.join(SCRIPT_DIR, '03_annotate_origin.py'),
           combined, annotated, sample_name,
           '--human_fasta', ref['human_peptide_fasta'],
           '--virus_fasta', ref['virus_peptide_fasta'],
           '--proteome_dir', ref.get('proteome_dir', '.')]
    if ref.get('proteome_fastas'):
        cmd += ['--proteome_files'] + list(ref['proteome_fastas'])
    run_step(cmd, log_file)
    outputs['annotated'] = annotated
    tracker.end_step("annotate origin")

    # ---- step 4: add R0 ----
    tracker.start_step(f"[{sample_name}] add R0 counts")
    with_r0 = os.path.join(results_dir, f"{sample_name}_R1R2_annotated_with_R0.csv")
    run_step([sys.executable, os.path.join(SCRIPT_DIR, '04_add_round0.py'),
              annotated, with_r0,
              '--r0_virus', r0['virus_counts'], '--r0_human', r0['human_counts']],
             log_file)
    outputs['with_r0'] = with_r0
    tracker.end_step("add R0 counts")

    # ---- step 5: filter matched ----
    tracker.start_step(f"[{sample_name}] filter matched peptides")
    matched = os.path.join(
        results_dir, f"{sample_name}_R1R2_annotated_with_R0_matched_filtered.csv")
    run_step([sys.executable, os.path.join(SCRIPT_DIR, '05_filter_matched.py'),
              with_r0, matched], log_file)
    outputs['matched_filtered'] = matched
    tracker.end_step("filter matched peptides")

    # ---- step 6: primers ----
    tracker.start_step(f"[{sample_name}] add primers")
    with_primer = os.path.join(
        results_dir,
        f"{sample_name}_R1R2_annotated_with_R0_matched_filtered_with_primer.csv")
    run_step([sys.executable, os.path.join(SCRIPT_DIR, '06_add_primers.py'),
              matched, with_primer, sample_name], log_file)
    outputs['with_primer'] = with_primer
    tracker.end_step("add primers")

    return outputs


def write_report(config_path, log_file, sample_outputs, total_time):
    first = next(iter(sample_outputs.values()))
    results_dir = os.path.dirname(first['with_primer'])
    report_path = os.path.join(results_dir, 'peptide_screen_processing_report.md')

    lines = [
        "# Peptide screen — processing report",
        "",
        f"- Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"- Config: `{config_path}`",
        f"- Total time: {timedelta(seconds=int(total_time))}",
        f"- Log: `{log_file}`",
        "",
        "## Per-sample outputs",
        "",
        "| Sample | Combined | Annotated | With R0 | Matched-filtered | Final (with primers) |",
        "|---|---|---|---|---|---|",
    ]
    for sample, outputs in sample_outputs.items():
        lines.append("| {s} | `{c}` | `{a}` | `{r}` | `{m}` | `{p}` |".format(
            s=sample,
            c=os.path.basename(outputs['combined']),
            a=os.path.basename(outputs['annotated']),
            r=os.path.basename(outputs['with_r0']),
            m=os.path.basename(outputs['matched_filtered']),
            p=os.path.basename(outputs['with_primer']),
        ))
    lines += [
        "",
        "## Pipeline steps",
        "",
        "1. `01_extract_peptide_counts.py` — FASTQ → `fullresult_*.csv` (nt, count, aa, length)",
        "2. `02_combine_rounds.py` — merge R1+R2 (+ percentages, r2_r1_ratio)",
        "3. `03_annotate_origin.py` — match peptide DNA to human/virus libraries, "
        "add fasta_name/origin/Gene",
        "4. `04_add_round0.py` — add count_round0 / count_pct_round0 / r1_r0_ratio",
        "5. `05_filter_matched.py` — keep origin in {human, virus}",
        "6. `06_add_primers.py` — sort by count_round2 desc, renumber, add primer-F/R",
        "",
        "## Final file columns",
        "",
        "`nt, aa, length, count_pct_round0, count_round0, r1_r0_ratio, count_pct_round1, "
        "count_round1, count_pct_round2, count_round2, r2_r1_ratio, fasta_name, origin, Gene, "
        "peptide_number, primer_F_name, primer-F, primer_R_name, primer-R`",
        "",
    ]
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f"\nReport written: {report_path}")
    return report_path


def main():
    ap = argparse.ArgumentParser(description="Run the Illumina peptide-only screen pipeline")
    ap.add_argument('--config', required=True,
                    help='YAML config (see configs/peptide_screen.yml.example)')
    args = ap.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    results_dir = cfg.get('output_dir', 'results')
    logs_dir = cfg.get('logs_dir', 'logs')
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(logs_dir, f"peptide_screen_{stamp}.log")

    samples = cfg['samples']
    steps_per_sample = 6
    total_steps = steps_per_sample * len(samples)
    tracker = ProgressTracker(total_steps, log_file)

    with open(log_file, 'w') as lf:
        lf.write("B27_single_TCR_screen — peptide screen\n")
        lf.write(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        lf.write(f"Config: {os.path.abspath(args.config)}\n")
        lf.write(f"Samples: {', '.join(samples.keys())}\n")
        lf.write("=" * 60 + "\n")

    start_time = time.time()
    sample_outputs = {}
    try:
        for sample_name, sample_cfg in samples.items():
            sample_outputs[sample_name] = run_sample(
                sample_name, sample_cfg, cfg, log_file, tracker)
    except Exception as e:
        msg = f"\nPIPELINE FAILED: {e}\nLog: {log_file}\n"
        print(msg)
        with open(log_file, 'a') as lf:
            lf.write(msg)
        sys.exit(1)

    total_time = time.time() - start_time
    report = write_report(args.config, log_file, sample_outputs, total_time)

    final = (f"\n{'='*60}\nAll samples processed!\n"
             f"Total time: {timedelta(seconds=int(total_time))}\n"
             f"Log: {log_file}\nReport: {report}\n{'='*60}")
    print(final)
    with open(log_file, 'a') as lf:
        lf.write(final)


if __name__ == "__main__":
    main()
