# -*- coding: utf-8 -*-
import os
import re
import logging
from collections import Counter
import tempfile
from utils import (read_fasta, reverse_complement, run_rnahybrid, read_strand_info,
                   calculate_gc_content, predict_structure_and_energy, Hairpin_structure,
                   parse_gene_info, energy_file_lock, global_print_lock,
                   _safe_up_left, fetch_transcript_region)
from config import XI_LIST, THETA_LIST, EXTEND_UP, EXTEND_DOWN, STRONG_SD_THRESHOLD, BASIC_SD_THRESHOLD, UTNI_ENERGY_FILE, RRNA_MOTIF_FOLDER, SD_MOTIF_SUFFIX

#UTNI_start_codon_regions = {}

def process_UTNI_gene(UTNI_gene, chromosome, strand_dict, genome_seq, species_name, energy_file_path, start_codon_pos_dict, overlap_gene_end):
    start_codon_pos = start_codon_pos_dict.get(UTNI_gene)
    if start_codon_pos is None:
        print(f"Warning: Start codon position not found for gene {UTNI_gene}")
        return False, None, 0, {}
    strand = strand_dict.get(UTNI_gene)
    if strand not in ('+', '-'):
        print(f"Warning: Unknown strand {strand} for {UTNI_gene}")
        return False, None, 0, {}

    safe_left = _safe_up_left(start_codon_pos, overlap_gene_end, strand)
    if safe_left is None:
        print(f"Warning: safe_left is None for {UTNI_gene}")
        return False, None, 0, {}

    raw_chr_seq = genome_seq[chromosome]
    full_RNA_binding_region, RNA_binding_region = fetch_transcript_region(
        raw_chr_seq=raw_chr_seq,
        start_codon_pos=start_codon_pos,
        safe_left=safe_left,
        strand=strand,
        gene_id=UTNI_gene,
        debug=False
    )
    if not full_RNA_binding_region:
        return False, None, 0, {}

    try:
        rRNA_end_sequence_file = os.path.join(RRNA_MOTIF_FOLDER, f"{species_name}{SD_MOTIF_SUFFIX}")
        RNA_binding_region = RNA_binding_region.replace('T', 'U')

        interaction_free_energy = None
        if len(RNA_binding_region.strip()) > 0:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_fa:
                tmp_fa.write(f">{UTNI_gene}\n{RNA_binding_region}\n")
                tmp_fa_path = tmp_fa.name
            interaction_free_energy = run_rnahybrid(tmp_fa_path, rRNA_end_sequence_file, XI_LIST, THETA_LIST)
            os.unlink(tmp_fa_path)

        UTNI_rbs_strength = "non_SD match"
        if interaction_free_energy is not None:
            _, _, energy = interaction_free_energy
            if energy is not None:
                if energy <= STRONG_SD_THRESHOLD:
                    UTNI_rbs_strength = "Strong SD"
                elif energy <= BASIC_SD_THRESHOLD:
                    UTNI_rbs_strength = "Basic SD"
                else:
                    UTNI_rbs_strength = "non_TC translation"

        with energy_file_lock:
            with open(energy_file_path, 'a') as ef:
                if interaction_free_energy is not None:
                    ef.write(f"{UTNI_gene}\t{interaction_free_energy}\t{UTNI_rbs_strength}\n")
                else:
                    ef.write(f"{UTNI_gene}\tNA\tnon_SD match\n")

        gene_result_dict = {
            'full_RNA_binding_region': (
                full_RNA_binding_region,          
                full_RNA_binding_region.replace('T', 'U'), 
                interaction_free_energy,        
                UTNI_rbs_strength                 
            ),
            'RNA_binding_region': RNA_binding_region  
        }
        full_rna_seq = full_RNA_binding_region.replace('T', 'U')
        structure, _ = predict_structure_and_energy(full_rna_seq)
#       print(f"Gene {UTNI_gene} RNA structure: {structure}")
        gc_content = calculate_gc_content(full_RNA_binding_region)
        if structure and Hairpin_structure(structure):
            return True, UTNI_rbs_strength, gc_content, gene_result_dict
        else:
            return False, UTNI_rbs_strength, gc_content, gene_result_dict
    except Exception as e:
        print(f"An error occurred in process_UTNI_gene {UTNI_gene}: {e}")
        return False, None, 0, {}

def count_UTNI_genes(gene_annotation_file, species_name, rRNA_folder, gff3_folder, fasta_folder):
    UTNI_gene_count = 0
    processed_gene_pairs = set()
    processed_genes = set()
    UTNI_genes = []
    UTNI_rbs_strength_counter = Counter()
    UTNI_gc_total = 0
    seq_dict = {}  
    gff3_filename = os.path.basename(gene_annotation_file)
    gff3_file = os.path.join(gff3_folder, gff3_filename)
    strand_dict = read_strand_info(gff3_file)
    gene_info_dict, start_codon_pos_dict = parse_gene_info(gff3_file)
    fasta_filename = os.path.basename(gene_annotation_file).replace('.gff3', '.fasta')
    genome_fasta = os.path.join(fasta_folder, fasta_filename)
    genome_seq = read_fasta(genome_fasta)
    interaction_energy = UTNI_ENERGY_FILE

    try:
        sorted_genes = sorted(gene_info_dict.items(), key=lambda item: (item[1][3], item[1][0]))
    except Exception as e:
        logging.error(f"[{species_name}] Sort gene list failed: {e}")
        return UTNI_gene_count, processed_gene_pairs, UTNI_genes, UTNI_rbs_strength_counter, seq_dict

    for i in range(len(sorted_genes) - 1):
        gene1, (start1, end1, strand1, chromosome1) = sorted_genes[i]
        gene2, (start2, end2, strand2, chromosome2) = sorted_genes[i + 1]

        if chromosome1 != chromosome2 or strand1 != strand2:
            continue
        distance = start2 - end1
        if 0 < distance < 20:
            gene_pair = tuple(sorted([gene1, gene2]))
            processed_gene_pairs.add(gene_pair)

            if strand1 == '+':
                UTNI_gene = gene2
                prev_stop = end1
            else:
                UTNI_gene = gene1
                prev_stop = end2
            if UTNI_gene in processed_genes:
                continue
    
            try:
                res, rbs_strength, gc_content, gene_result_dict = process_UTNI_gene(
                    UTNI_gene=UTNI_gene,
                    chromosome=chromosome1,
                    strand_dict=strand_dict,
                    genome_seq=genome_seq,
                    species_name=species_name,
                    energy_file_path=interaction_energy,
                    start_codon_pos_dict=start_codon_pos_dict,
                    overlap_gene_end=prev_stop
                )
                if res:
                    UTNI_gene_count += 1
                    processed_genes.add(UTNI_gene)
                    UTNI_genes.append(UTNI_gene)
                    seq_dict[UTNI_gene] = gene_result_dict
                    UTNI_gc_total += gc_content
                    UTNI_rbs_strength_counter[rbs_strength] += 1
            except Exception as e:
                logging.error(f"[{species_name}] Gene {UTNI_gene} processing failed, skip. Error: {e}")
                continue
    
    return UTNI_gene_count, processed_gene_pairs, UTNI_genes, UTNI_rbs_strength_counter, seq_dict