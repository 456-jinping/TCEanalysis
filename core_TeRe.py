# -*- coding: utf-8 -*-
import os
import re
import logging
from collections import Counter
import tempfile
from utils import (read_fasta, reverse_complement, run_rnahybrid, read_strand_info,
                   energy_file_lock, global_print_lock, _safe_up_left, fetch_transcript_region)
from config import (XI_LIST, THETA_LIST, EXTEND_UP, EXTEND_DOWN, STRONG_SD_THRESHOLD,
                    BASIC_SD_THRESHOLD, TERE_ENERGY_FILE, RRNA_MOTIF_FOLDER, SD_MOTIF_SUFFIX)

gene_id_pattern = re.compile(r'ID=([^;]+)')
gene_name_pattern = re.compile(r'gene=([^;]+)')

def parse_gff(gene_annotation_file):
    strand_dict = {}
    start_codon_pos_dict = {}
    gene_end_dict = {}  
    gene_start_dict = {}
    with open(gene_annotation_file, 'r', encoding="utf-8") as ann_file:
        for line in ann_file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            attributes = cols[8].split(';')
            gene_id = None
            chromosome = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]
            for attr in attributes:
                attr = attr.strip()
                if attr.startswith('ID='):
                    gene_id = attr.split('=')[1].strip('"')
            if gene_id:
                strand_dict[gene_id] = strand
                gene_start_dict[gene_id] = start
                gene_end_dict[gene_id] = end
                if strand == '-':
                    start_codon_pos_dict[gene_id] = end
                else:
                    start_codon_pos_dict[gene_id] = start
    return strand_dict, start_codon_pos_dict, gene_end_dict, gene_start_dict

def scan_bed_full(bed_file):
    bed_lines = []
    bed_gene_map = dict()
    raw_pairs = []
    with open(bed_file, 'r', encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 8:
                continue
            g1_info = cols[3]
            g2_info = cols[7]
            m1 = gene_id_pattern.search(g1_info) or gene_name_pattern.search(g1_info)
            m2 = gene_id_pattern.search(g2_info) or gene_name_pattern.search(g2_info)
            if not m1 or not m2:
                continue
            g1 = m1.group(1).strip('"')
            g2 = m2.group(1).strip('"')
            if g1 == g2:
                continue
            bed_gene_map[g1] = line
            bed_gene_map[g2] = line
            raw_pairs.append((g1, g2))
            bed_lines.append(line)
    return bed_lines, bed_gene_map, raw_pairs

def process_TeRe_gene(TeRe_gene, bed_line, genome_seq, strand_dict, start_codon_pos_dict, species_name, chr_comp_cache, overlap_gene_end):
    start_codon_pos = start_codon_pos_dict.get(TeRe_gene)
    if start_codon_pos is None:
        return "", ""
    chromosome = bed_line.split('\t')[0]
    if chromosome not in genome_seq:
        return "", ""
    raw_chr_seq = genome_seq[chromosome]
    strand = strand_dict[TeRe_gene]
    safe_left = _safe_up_left(start_codon_pos, overlap_gene_end, strand)
    full_RNA_binding_region, RNA_binding_region = fetch_transcript_region(raw_chr_seq, start_codon_pos, safe_left, strand, gene_id=TeRe_gene, debug=False)
    return full_RNA_binding_region, RNA_binding_region

def count_overlapping_genes(bed_file, genome_fasta, gene_annotation_file, species_name):
    strand_dict, start_codon_pos_dict, gene_end_dict, gene_start_dict = parse_gff(gene_annotation_file)
    genome_seq = read_fasta(genome_fasta)
    bed_lines, bed_gene_map, raw_gene_pairs = scan_bed_full(bed_file)
    rRNA_motif_path = os.path.join(RRNA_MOTIF_FOLDER, f"{species_name}{SD_MOTIF_SUFFIX}")

    processed_gene_pairs = set()
    start_codon_regions = {}
    overlap_rbs_strength_counter = Counter()
    chr_comp_cache = dict()
    gene_pair_count = 0
    overlap_genes = set()  

    for g1, g2 in raw_gene_pairs:
        s1 = strand_dict.get(g1)
        s2 = strand_dict.get(g2)
        if not ((s1 == '+' and s2 == '+') or (s1 == '-' and s2 == '-')):
            continue
    
        if g1 not in gene_start_dict or g2 not in gene_start_dict:
            continue
        g1s = gene_start_dict[g1]
        g1e = gene_end_dict[g1]
        g2s = gene_start_dict[g2]
        g2e = gene_end_dict[g2]
    
        overlap_start = max(g1s, g2s)
        overlap_end = min(g1e, g2e)
        overlap_len = overlap_end - overlap_start + 1
        if overlap_len <= 0:
            continue
    
        pair_key = tuple(sorted([g1, g2]))
        if pair_key in processed_gene_pairs:
            continue
        processed_gene_pairs.add(pair_key)
        gene_pair_count += 1
#       overlap_genes.add(g1)
#       overlap_genes.add(g2)

        s = s1
        for TeRe_gene in [g1, g2]:
            if s == '+' and TeRe_gene == g1:
                continue
            if s == '-' and TeRe_gene == g2:
                continue
            if TeRe_gene in start_codon_regions:
                continue
            bed_line = bed_gene_map[TeRe_gene]
            if TeRe_gene == g1:
                overlap_gene_end = gene_end_dict.get(g2, -1)
            else:
                overlap_gene_end = gene_end_dict.get(g1, -1)
            full_seq, sd_seq = process_TeRe_gene(TeRe_gene, bed_line, genome_seq, strand_dict, start_codon_pos_dict, species_name, chr_comp_cache, overlap_gene_end)
            if len(full_seq.strip()) == 0:
                continue
            interaction_free_energy = None
            if len(sd_seq.strip()) > 0:
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmpfa:
                    tmpfa.write(f">{TeRe_gene}\n{sd_seq.replace('T', 'U')}\n")
                    tmp_path = tmpfa.name
                interaction_free_energy = run_rnahybrid(tmp_path, rRNA_motif_path, XI_LIST, THETA_LIST)
                os.unlink(tmp_path)

            rbs_strength = "non_SD match"
            if interaction_free_energy is not None:
                _, _, energy = interaction_free_energy
                if energy is not None:
                    if energy <= STRONG_SD_THRESHOLD:
                        rbs_strength = "Strong SD"
                    elif energy <= BASIC_SD_THRESHOLD:
                        rbs_strength = "Basic SD"
                    else:
                        rbs_strength = "non_TC translation"
            overlap_rbs_strength_counter[rbs_strength] += 1
            with energy_file_lock:
                with open(TERE_ENERGY_FILE, 'a', encoding="utf-8") as ef:
                    if interaction_free_energy is not None:
                        ef.write(f"{TeRe_gene}\t{interaction_free_energy}\t{rbs_strength}\n")
                    else:
                        ef.write(f"{TeRe_gene}\tNA\tnon_SD match\n")
            start_codon_regions[TeRe_gene] = {
                'full_RNA_binding_region': (
                    full_seq,
                    full_seq.replace('T', 'U'),
                    interaction_free_energy,
                    rbs_strength
                ),
                'RNA_binding_region': sd_seq,
            }
            
    return start_codon_regions, processed_gene_pairs, gene_pair_count, overlap_rbs_strength_counter