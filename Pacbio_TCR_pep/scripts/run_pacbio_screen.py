#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Driver (pacbio_tcr_peptide): run the full PacBio TCR+peptide screen from config:

    Stage 1  per round: BAM -> FASTQ (samtools bam2fq), then
             01_extract_tcr_peptide.py (minimap2 TCR alignment + peptide extraction)
    Stage 2  per round: 02_final_qc.py -> final_QC.csv
    Stage 3  03_solo_analysis.py on the two final_QC tables ->
             filtered_data_clean_with_primer.csv

Every command and its output is appended to a run log; a Markdown processing
report is written at the end. Progress and estimated remaining time are shown
after each step.

Usage:
    python run_pacbio_screen.py --config configs/my_pacbio.yml
"""

import os
import sys
import glob
import time
import argparse
import subprocess
from datetime import datetime, timedelta

import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def fmt_time(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    if seconds < 3600:
        return f"{seconds/60:.1f}min"
    return f"{seconds/3600:.2f}h"


def log(msg, log_file=None):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line)
    if log_file:
        with open(log_file, 'a') as f:
            f.write(line + '\n')


def run_cmd(cmd, log_file):
    printable = ' '.join(str(c) for c in cmd)
    log(f"[CMD] {printable}", log_file)
    start = time.time()
    proc = subprocess.Popen([str(c) for c in cmd],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in proc.stdout:
        print(line, end='')
        with open(log_file, 'a') as f:
            f.write(line)
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed (exit {proc.returncode}): {printable}")
    log(f"[DONE] {printable} ({fmt_time(time.time()-start)})", log_file)


class ProgressTracker:
    def __init__(self, total_steps, log_file):
        self.total_steps = total_steps
        self.current_step = 0
        self.start_time = time.time()
        self.log_file = log_file

    def start_step(self, name):
        self.current_step += 1
        self.step_start = time.time()
        log(f"\n{'='*70}\nStep {self.current_step}/{self.total_steps}: {name}\n{'='*70}",
            self.log_file)

    def end_step(self, name):
        elapsed = time.time() - self.step_start
        total = time.time() - self.start_time
        avg = total / self.current_step
        remaining = avg * (self.total_steps - self.current_step)
        log(f"Step done: {name} | this step {fmt_time(elapsed)} | "
            f"total {fmt_time(total)} | estimated remaining {fmt_time(remaining)}\n{'='*70}",
            self.log_file)


def bam_to_fastq(bam_file, output_fastq, log_file):
    log(f"  BAM -> FASTQ: {bam_file}", log_file)
    cmd = f"samtools bam2fq -0 {output_fastq} {bam_file}"
    subprocess.run(cmd, shell=True, check=True)


def process_round(round_cfg, cfg, log_file, data_dir):
    """Stage 1 for one round: BAM -> FASTQ -> extraction. Returns result prefix."""
    name = round_cfg['name']
    round_dir = os.path.join(data_dir, name)
    output_dir = os.path.join(round_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)

    fastq_file = os.path.join(output_dir, f'{name}.hifi_reads.fastq')
    result_prefix = os.path.join(output_dir, name)
    result_csv = f"{result_prefix}_result.csv"

    if not os.path.exists(fastq_file):
        if round_cfg.get('fastq'):
            log(f"  Using provided FASTQ: {round_cfg['fastq']}", log_file)
            fastq_file = round_cfg['fastq']
        else:
            bam_files = glob.glob(round_cfg['bam_glob'])
            if not bam_files:
                raise FileNotFoundError(
                    f"No BAM files found for round {name}: {round_cfg['bam_glob']}")
            bam_to_fastq(bam_files[0], fastq_file, log_file)
    else:
        log(f"  FASTQ already exists, skipping conversion: {fastq_file}", log_file)

    if not os.path.exists(result_csv):
        extract = cfg.get('extract', {})
        run_cmd([sys.executable, os.path.join(SCRIPT_DIR, '01_extract_tcr_peptide.py'),
                 fastq_file, cfg['reference']['tcr_reference'], result_prefix,
                 '--threads', extract.get('threads', 8),
                 '--min_len', extract.get('min_insert_len', 1800),
                 '--max_len', extract.get('max_insert_len', 2500)], log_file)
    else:
        log(f"  result.csv already exists, skipping extraction: {result_csv}", log_file)

    return result_prefix, fastq_file


def main():
    ap = argparse.ArgumentParser(description="Run the PacBio TCR+peptide screen pipeline")
    ap.add_argument('--config', required=True,
                    help='YAML config (see configs/pacbio_screen.yml.example)')
    args = ap.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    sample_name = cfg['sample_name']
    analysis_dir = cfg.get('analysis_dir', 'analysis')
    output_dir = cfg.get('output_dir', 'output')
    logs_dir = cfg.get('logs_dir', 'logs')
    for d in (analysis_dir, output_dir, logs_dir):
        os.makedirs(d, exist_ok=True)

    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(logs_dir, f"pacbio_screen_{stamp}.log")

    rounds = cfg['rounds']
    total_steps = len(rounds) * 2 + 1   # stage1+stage2 per round, + solo analysis
    tracker = ProgressTracker(total_steps, log_file)

    log("=" * 70, log_file)
    log("B27_single_TCR_screen — PacBio TCR+peptide screen", log_file)
    log(f"Config: {os.path.abspath(args.config)}", log_file)
    log(f"Sample: {sample_name} | Rounds: {[r['name'] for r in rounds]}", log_file)
    log("=" * 70, log_file)

    start_time = time.time()
    result_prefixes, final_qc_files = [], []

    try:
        # ---- Stage 1 ----
        for i, rnd in enumerate(rounds, 1):
            tracker.start_step(f"Stage 1 [{i}/{len(rounds)}]: {rnd['name']} "
                               f"(BAM->FASTQ + TCR/peptide extraction)")
            t0 = time.time()
            prefix, fastq_file = process_round(rnd, cfg, log_file, output_dir)
            bam_file = (glob.glob(rnd['bam_glob'])[0]
                        if rnd.get('bam_glob') and glob.glob(rnd['bam_glob']) else None)
            log(f"  BAM size: {os.path.getsize(bam_file)/1024**3:.2f} GB"
                if bam_file else "  (no BAM size info)", log_file)
            log(f"  FASTQ size: {os.path.getsize(fastq_file)/1024**3:.2f} GB"
                if os.path.exists(fastq_file) else "", log_file)
            result_prefixes.append(prefix)
            tracker.end_step(f"Stage 1 [{i}/{len(rounds)}]: {rnd['name']}")

        # ---- Stage 2 ----
        for i, prefix in enumerate(result_prefixes, 1):
            name = os.path.basename(prefix)
            tracker.start_step(f"Stage 2 [{i}/{len(result_prefixes)}]: Final QC {name}")
            result_csv = f"{prefix}_result.csv"
            output_dir_round = os.path.join(output_dir, name, 'output')
            final_qc = run_qc_step(result_csv, output_dir_round, name, log_file)
            final_qc_files.append(final_qc)
            tracker.end_step(f"Stage 2 [{i}/{len(result_prefixes)}]: Final QC {name}")

        # ---- Stage 3 ----
        tracker.start_step("Stage 3: solo two-round analysis + primers")
        if len(final_qc_files) < 2:
            raise RuntimeError("Need final_QC.csv from at least 2 rounds for solo analysis")
        analysis_cfg = cfg.get('analysis', {})
        run_cmd([sys.executable, os.path.join(SCRIPT_DIR, '03_solo_analysis.py'),
                 final_qc_files[0], final_qc_files[1],
                 '--output_dir', analysis_dir,
                 '--sample_name', sample_name,
                 '--human_proteome', cfg['reference']['human_proteome'],
                 '--virus_proteome', cfg['reference']['virus_proteome'],
                 '--allowed_tcr_numbers', analysis_cfg.get('allowed_tcr_numbers', '1-49,101'),
                 '--fold_change_min', analysis_cfg.get('fold_change_min', 0.1),
                 '--fold_change_max', analysis_cfg.get('fold_change_max', 60),
                 '--min_count_ratio', analysis_cfg.get('min_count_ratio', 0.05),
                 '--min_peptide_len', analysis_cfg.get('min_peptide_len', 8),
                 '--special_peptides'] + list(analysis_cfg.get('special_peptides', [])) +
                ['--remove_keywords'] + list(analysis_cfg.get('remove_peptide_keywords', [])),
                log_file)
        tracker.end_step("Stage 3: solo two-round analysis + primers")

    except Exception as e:
        log(f"\nPIPELINE FAILED: {e}\nLog: {log_file}", log_file)
        sys.exit(1)

    total_time = time.time() - start_time
    primer_file = os.path.join(analysis_dir, 'filtered_data_clean_with_primer.csv')

    # ---- Markdown report ----
    report_path = os.path.join(output_dir, 'pacbio_screen_processing_report.md')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(f"# PacBio TCR+peptide screen — processing report\n\n")
        f.write(f"- Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"- Config: `{os.path.abspath(args.config)}`\n")
        f.write(f"- Sample: {sample_name}\n")
        f.write(f"- Total time: {timedelta(seconds=int(total_time))}\n")
        f.write(f"- Log: `{os.path.abspath(log_file)}`\n\n")
        f.write(f"## Pipeline\n\n")
        f.write(f"```text\nBAM (per round)\n  -> samtools bam2fq\n")
        f.write(f"  -> 01_extract_tcr_peptide.py (regex extraction + minimap2 TCR alignment)\n")
        f.write(f"  -> <round>_result.csv / _result_qc.csv\n")
        f.write(f"  -> 02_final_qc.py -> final_QC.csv (per round)\n")
        f.write(f"  -> 03_solo_analysis.py (merge + origin + fold change + dedup + clean filter)\n")
        f.write(f"  -> analysis/filtered_data_clean_with_primer.csv\n```\n\n")
        f.write(f"## Outputs\n\n")
        f.write(f"| Priority | File | Use |\n|---|---|---|\n")
        f.write(f"| 1 | `analysis/filtered_data_clean_with_primer.csv` | oligo synthesis |\n")
        f.write(f"| 2 | `analysis/filtered_data_clean.csv` | analysis / plotting |\n")
        f.write(f"| 3 | `analysis/{sample_name}_analysis_result_solo_filtered.csv` | deduped full data |\n")
        f.write(f"| 4 | `analysis/{sample_name}_analysis_result_with_origin.csv` | full data with origin |\n\n")
        f.write(f"## Key parameters\n\n")
        ac = cfg.get('analysis', {})
        f.write(f"- fold change window: [{ac.get('fold_change_min', 0.1)}, "
                f"{ac.get('fold_change_max', 60)}]\n")
        f.write(f"- per-peptide min count ratio: {ac.get('min_count_ratio', 0.05)}\n")
        f.write(f"- allowed TCR numbers: {ac.get('allowed_tcr_numbers', '1-49,101')}\n")
        f.write(f"- removed peptide keywords: {ac.get('remove_peptide_keywords', [])}\n")

    log("\n" + "=" * 70, log_file)
    log("All done!", log_file)
    log(f"Total time: {fmt_time(total_time)}", log_file)
    log(f"Final result: {os.path.abspath(primer_file)}", log_file)
    log(f"Report: {os.path.abspath(report_path)}", log_file)
    log(f"Log: {os.path.abspath(log_file)}", log_file)
    log("=" * 70, log_file)


def run_qc_step(result_csv, output_dir, name, log_file):
    run_cmd([sys.executable, os.path.join(SCRIPT_DIR, '02_final_qc.py'),
             result_csv, output_dir, name], log_file)
    return os.path.join(output_dir, 'final_QC.csv')


if __name__ == "__main__":
    main()
