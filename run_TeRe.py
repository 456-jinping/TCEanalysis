# -*- coding: utf-8 -*-
import logging
import os
import csv
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from core_TeRe import count_overlapping_genes
from config import FASTA_FOLDER, GFF3_FOLDER, BED_FOLDER, TERE_FASTA_OUT, TERE_MAX_WORKERS, ERROR_REPORT, TERE_STAT_EXCEL, BED_SUFFIX, SD_MOTIF_SUFFIX, RRNA_MOTIF_FOLDER, TERE_DONE_RECORD
from utils import safe_div, read_fasta, count_all_genes, global_print_lock, load_finished_species, mark_species_done, calculate_gc_content
import threading
csv_write_lock = threading.Lock()
energy_write_lock = threading.Lock()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
results = []

TERE_STAT_CSV = TERE_STAT_EXCEL.replace('.xlsx', '.csv')

def append_single_row_to_csv(out_csv, row_dict):
    columns = ['Species', 'All Gene Count', 'TeRe TCEs Count', 'TeRe TCEs Ratio',
               'TeRe TCEs with SD', 'TeRe TCEs with SD%', 'Strong SD Count (Overlap)',
               'TeRe TCEs with str_SD%', 'Basic SD Count (Overlap)']
    with csv_write_lock:
        file_exists = os.path.isfile(out_csv)
        with open(out_csv, 'a', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=columns)
            if not file_exists:
                writer.writeheader()
            writer.writerow(row_dict)
            
def process_file(filename):
    global results
    if not filename.endswith('.bed'):
        return filename
    bed_file = os.path.join(BED_FOLDER, filename)
    species_name = filename.replace(BED_SUFFIX, "")
    genome_fasta = os.path.join(FASTA_FOLDER, f'{species_name}.fasta')
    gene_annotation_file = os.path.join(GFF3_FOLDER, f'{species_name}.gff3')
    rRNA_end_sequence_file = os.path.join(RRNA_MOTIF_FOLDER, f"{species_name}{SD_MOTIF_SUFFIX}")
    req_files = [bed_file, genome_fasta, gene_annotation_file, rRNA_end_sequence_file]
    missing = [f for f in req_files if not os.path.exists(f)]
    if missing:
        logging.error(f"[{species_name}] Missing files: {missing}")
        return filename
    try:
        if not os.path.exists(TERE_FASTA_OUT):
            os.makedirs(TERE_FASTA_OUT)
#       sp_fasta_path = os.path.join(TERE_FASTA_OUT, f"{species_name}_overlap_sequences.fasta")
        sd_fasta_dict = read_fasta(rRNA_end_sequence_file)
        if not sd_fasta_dict:
            logging.error(f"[{species_name}] SD motif fasta empty")
            return filename
        rRNA_end_sequence = list(sd_fasta_dict.values())[0]
        all_gene_count = count_all_genes(gene_annotation_file)
        start_codon_regions, new_processed_gene_pairs, gene_pair_count, overlap_rbs_strength_counter = count_overlapping_genes(
            bed_file, genome_fasta, gene_annotation_file, species_name
        )
        overlap_gene_count = len(new_processed_gene_pairs)
        overlap_ratio = safe_div(overlap_gene_count, all_gene_count)
        strong_sd_count = overlap_rbs_strength_counter.get("Strong SD", 0)
        basic_sd_count = overlap_rbs_strength_counter.get("Basic SD", 0)
        Overlap_Gene_with_SD = strong_sd_count + basic_sd_count
        Overlap_Gene_with_SD_Ratio = safe_div(Overlap_Gene_with_SD, overlap_gene_count)
        Overlap_Gene_with_strSD_Ratio = safe_div(strong_sd_count, overlap_gene_count)
        result = {
            'Species': species_name,
            'All Gene Count': all_gene_count,
            'TeRe TCEs Count': overlap_gene_count,
            'TeRe TCEs Ratio': overlap_ratio,
            'TeRe TCEs with SD': Overlap_Gene_with_SD,
            'TeRe TCEs with SD%': Overlap_Gene_with_SD_Ratio,
            'Strong SD Count (Overlap)': strong_sd_count,
            'TeRe TCEs with str_SD%': Overlap_Gene_with_strSD_Ratio,
            'Basic SD Count (Overlap)': basic_sd_count
        }
        with global_print_lock:
            results.append(result)
            append_single_row_to_csv(TERE_STAT_CSV, result)
            
        if start_codon_regions:
            sp_fasta_path = os.path.join(TERE_FASTA_OUT, f"{species_name}.fasta")
            sp_meta_path = os.path.join(TERE_FASTA_OUT, f"{species_name}.txt")
            with open(sp_fasta_path, "w", encoding="utf-8") as fa_out, open(sp_meta_path, "w", encoding="utf-8") as meta_out:
                meta_out.write("gene_id\tfull_dna_seq\tgc_ratio\txi\ttheta\tmfe_energy\trbs_strength\n")
                for gid, data in start_codon_regions.items():
                    dna_seq = data["full_RNA_binding_region"][0]
                    fa_out.write(f">{gid}\n{dna_seq}\n")
                    full_u, xi, theta, mfe = data["full_RNA_binding_region"]
                    rbs_strength = data["full_RNA_binding_region"][3]
                    gc_ratio = calculate_gc_content(dna_seq) / 100
                    if mfe is None:
                        meta_out.write(f"{gid}\t{dna_seq}\t{gc_ratio:.4f}\tNA\tNA\tNA\t{rbs_strength}\n")
                    else:
                        meta_out.write(f"{gid}\t{dna_seq}\t{gc_ratio:.4f}\t{xi}\t{theta}\t{mfe}\t{rbs_strength}\n")
            logging.info(f"[{species_name}] FASTA + meta txt saved to {TERE_FASTA_OUT}")
        else:
            logging.warning(f"[{species_name}] No overlap TeRe genes found, skip fasta & meta txt")
            
        logging.info(f"[{species_name}] Process finished, TeRe count={overlap_gene_count}")
        mark_species_done(TERE_DONE_RECORD, species_name)
        return None
    except Exception as e:
        logging.error(f"[{species_name}] Process failed: {str(e)}", exc_info=True)
        return filename

if __name__ == '__main__':
    results.clear()
    done_set = load_finished_species(TERE_DONE_RECORD)
    all_bed_files = [f for f in os.listdir(BED_FOLDER) if f.endswith(".bed")]

    todo_bed = []
    for bed_name in all_bed_files:
        sp = bed_name.replace(BED_SUFFIX, "")
        if sp not in done_set:
            todo_bed.append(bed_name)

    print(f"Bed_file：{len(all_bed_files)} | done_species：{len(done_set)} | need_to_process：{len(todo_bed)}")
    if len(todo_bed) == 0:
        print("All genomes have processed, exit！")
        exit()

    unprocessed_species = set()
    with ThreadPoolExecutor(max_workers=TERE_MAX_WORKERS) as executor:
        futures = [executor.submit(process_file, fn) for fn in todo_bed]
        for fut in futures:
            err_file = fut.result()
            if err_file is not None:
                unprocessed_species.add(err_file)
    if unprocessed_species:
        with open(ERROR_REPORT, 'w', encoding="utf-8") as f:
            f.write("Unprocessed bed files / missing resources:\n")
            for fname in sorted(unprocessed_species):
                f.write(f"{fname}\n")
    if results:
        df = pd.DataFrame(results)
        df.to_excel(TERE_STAT_EXCEL, index=False)
        print(f"Final Excel saved to {TERE_STAT_EXCEL}")
    else:
        print("No results to save, Excel not generated.")
        
    df = pd.DataFrame(results)
    df.to_excel(TERE_STAT_EXCEL, index=False)
    print("All overlap gene analysis completed.")
    print(f"Statistics saved to {TERE_STAT_EXCEL}")
    print(f"TeRe gene fasta sequences saved in: {TERE_FASTA_OUT}")