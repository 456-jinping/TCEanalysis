# -*- coding: utf-8 -*-
import os
import tempfile
import logging
from collections import Counter
from utils import read_fasta, reverse_complement, run_rnahybrid, read_strand_info, parse_gene_info, safe_div
from config import XI_LIST, THETA_LIST, STRONG_SD_THRESHOLD, BASIC_SD_THRESHOLD, RRNA_MOTIF_FOLDER, LEADING_EXTEND_UP, LEADING_LONG_GAP, SD_MOTIF_SUFFIX

def process_leading_gene(leading_gene, chromosome, strand_dict, genome_seq, rev_comp_cache, species_name, energy_file, start_codon_pos_dict, leading_rbs_strength_counter):
    start_codon_pos = start_codon_pos_dict.get(leading_gene)
    if start_codon_pos is None:
        logging.warning(f"Start codon position not found for gene {leading_gene}")

        return {}

    RNA_binding_region = ""
    if strand_dict.get(leading_gene) == '+':
        region_start = max(0, start_codon_pos - LEADING_EXTEND_UP)
        region_end = start_codon_pos
        RNA_binding_region = genome_seq[chromosome][region_start:region_end]
    else:
        complementary_seq = rev_comp_cache[chromosome]
        complementary_start_codon_pos = len(complementary_seq) - start_codon_pos
        region_start = max(0, complementary_start_codon_pos - LEADING_EXTEND_UP)
        region_end = complementary_start_codon_pos
        RNA_binding_region = complementary_seq[region_start:region_end]

    rRNA_end_sequence_file = os.path.join(RRNA_MOTIF_FOLDER, f"{species_name}{SD_MOTIF_SUFFIX}")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmpfa:
        tmpfa.write(f">{leading_gene}\n{RNA_binding_region.replace('T', 'U')}\n")
        tmp_path = tmpfa.name
    interaction_free_energy = run_rnahybrid(tmp_path, rRNA_end_sequence_file, XI_LIST, THETA_LIST)
    os.unlink(tmp_path)

    leading_rbs_strength = "non_SD match"
    if interaction_free_energy is not None:
        _, _, energy = interaction_free_energy
        if energy <= STRONG_SD_THRESHOLD:
            leading_rbs_strength = "Strong SD"
        elif energy <= BASIC_SD_THRESHOLD:
            leading_rbs_strength = "Basic SD"
        else:
            leading_rbs_strength = "non_ribo translation"
    leading_rbs_strength_counter[leading_rbs_strength] += 1

    if interaction_free_energy is not None:
        energy_file.write(f"{leading_gene}\t{interaction_free_energy}\t{leading_rbs_strength}\n")
    else:
        energy_file.write(f"{leading_gene}\tNA\tnon_SD match\n")

    region_data = (
        RNA_binding_region,
        RNA_binding_region.replace('T', 'U'),
        interaction_free_energy,
        leading_rbs_strength
    )
    single_gene_dict = {}
    single_gene_dict[leading_gene] = {'RNA_binding_region': region_data}
    return single_gene_dict

def count_leading_genes(gene_annotation_file, species_name, gff3_folder, fasta_folder, species_energy_path):
    leading_gene_count = 0
    processed_gene_pairs = set()
    processed_genes = set()
    leading_rbs_strength_counter = Counter()

    gff3_filename = os.path.basename(gene_annotation_file)
    gff3_file = os.path.join(gff3_folder, gff3_filename)
    strand_dict = read_strand_info(gff3_file)
    gene_info_dict, start_codon_pos_dict = parse_gene_info(gene_annotation_file)

    fasta_filename = os.path.basename(gene_annotation_file).replace('.gff3', '.fasta')
    genome_fasta = os.path.join(fasta_folder, fasta_filename)
    genome_seq = read_fasta(genome_fasta)

    rev_comp_cache = {}
    for chrom, seq in genome_seq.items():
        rev_comp_cache[chrom] = reverse_complement(seq)

    all_leading_regions = {}
    with open(species_energy_path, 'a', encoding="utf-8") as energy_file:
        sorted_genes = sorted(gene_info_dict.items(), key=lambda item: (item[1][3], item[1][0]))
        def _proc(gene, chrom):
            nonlocal leading_gene_count
            if gene not in processed_genes:
                leading_gene_count += 1
                processed_genes.add(gene)
                rd = process_leading_gene(
                    gene, chrom, strand_dict, genome_seq, rev_comp_cache,
                    species_name, energy_file, start_codon_pos_dict, leading_rbs_strength_counter
                )
                all_leading_regions.update(rd)

        for i in range(len(sorted_genes) - 1):
            gene1, (start1, end1, strand1, chromosome1) = sorted_genes[i]
            gene2, (start2, end2, strand2, chromosome2) = sorted_genes[i + 1]
            if chromosome1 != chromosome2:
                continue

            if strand1 == '+' and strand2 == '+':
                gap = start2 - end1
                if gap > LEADING_LONG_GAP:
                    _proc(gene2, chromosome1)

            elif strand1 == '-' and strand2 == '-':
                gap = start2 - end1
                if gap > LEADING_LONG_GAP:
                    _proc(gene1, chromosome1)

            elif strand1 == '-' and strand2 == '+':
                _proc(gene1, chromosome1)
                _proc(gene2, chromosome2)

            gene_pair = tuple(sorted([gene1, gene2]))
            processed_gene_pairs.add(gene_pair)

    return leading_gene_count, processed_gene_pairs, all_leading_regions, leading_rbs_strength_counter
