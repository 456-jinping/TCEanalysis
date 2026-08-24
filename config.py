# -*- coding: utf-8 -*-
# ===================== Global shared paths =====================
FASTA_FOLDER = "/path/to/genenome/fasta/file/"
GFF3_FOLDER = "/path/to/genenome/gff3/file/"
RRNA_MOTIF_FOLDER = "/path/to//16s/rRNA/end/sequence/file/"

# RNAhybrid parameters
XI_LIST = [0.1, 0.2, 0.3, 0.4, 0.5]
THETA_LIST = [0.8, 0.9, 1.0, 1.1, 1.2]
# SD strength thresholds
STRONG_SD_THRESHOLD = -8.4
BASIC_SD_THRESHOLD = -3.5

EXTEND_UP = 30
EXTEND_DOWN = 15

SD_MOTIF_SUFFIX = "_sd_motifs.fasta"

# ===================== TeRe specific configuration =====================
BED_FOLDER = "/path/to/genenome/bed/file/" # use bedtools
TERE_FASTA_OUT = "./All_TeRe_genes/"
BED_SUFFIX = ".gene.overlap.bed"
TERE_ENERGY_FILE = "TeRe_energy.output"
TERE_DONE_RECORD = "TeRe_finished.txt"
TERE_STAT_EXCEL = "TeRe_TCEs.xlsx"
TERE_MAX_WORKERS = 20

# ===================== UTNI specific configuration =====================
UTNI_FASTA_OUT = "./All_UTNI_genes/"
UTNI_ENERGY_FILE = "UTNI_energy.output"
UTNI_DONE_RECORD = "UTNI_finished.txt"
UTNI_STAT_EXCEL = "UTNI_TCEs.xlsx"
UTNI_MAX_WORKERS = 20

# ===================== Leading specific configuration =====================
LEADING_FASTA_OUT = "./All_Leading_genes/"
LEADING_ENERGY_FILE = "Leadinggene_energy.output"
LEADING_DONE_RECORD = "Leadinggene_finished.txt"
LEADING_STAT_EXCEL = "Leading_genes.xlsx"
LEADING_MAX_WORKERS = 20
LEADING_EXTEND_UP = 30
LEADING_LONG_GAP = 200

# ===================== Global common output file =====================
ERROR_REPORT = "error_report.txt"

# Checkpoint file for resume‑from‑breakpoint
DONE_RECORD = "done_species.txt"