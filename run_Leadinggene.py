# -*- coding: utf-8 -*-
import os
import csv
import logging
import threading
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from utils import safe_div, read_fasta, count_all_genes, global_print_lock, load_finished_species, mark_species_done, calculate_gc_content
from core_Leadinggene import count_leading_genes
from config import FASTA_FOLDER, GFF3_FOLDER, RRNA_MOTIF_FOLDER, LEADING_FASTA_OUT, ERROR_REPORT, LEADING_STAT_EXCEL, SD_MOTIF_SUFFIX, LEADING_MAX_WORKERS, LEADING_DONE_RECORD

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
csv_write_lock = threading.Lock()
LEADING_STAT_CSV = LEADING_STAT_EXCEL.replace('.xlsx', '.csv')

def append_single_row_to_csv(out_csv, row_dict):
    columns = [
        'Species', 'All Gene Count', 'Leading Gene Count', 'Leading Gene Ratio',
        'Leading Gene with SD', 'Leading Gene with SD%', 'Strong SD Count (Leading)',
        'Leading Gene with str_SD%', 'Basic SD Count (Leading)'
    ]
    with csv_write_lock:
        file_exists = os.path.isfile(out_csv)
        with open(out_csv, 'a', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=columns)
            if not file_exists:
                writer.writeheader()
            writer.writerow(row_dict)

def process_file(filename):
    if not filename.endswith('.gff3'):
        return filename
    species_name = filename.rsplit('.', 1)[0]
    genome_fasta = os.path.join(FASTA_FOLDER, f'{species_name}.fasta')
    gene_annotation_file = os.path.join(GFF3_FOLDER, f'{species_name}.gff3')
    rRNA_end_sequence_file = os.path.join(RRNA_MOTIF_FOLDER, f"{species_name}{SD_MOTIF_SUFFIX}")
    req_files = [genome_fasta, gene_annotation_file, rRNA_end_sequence_file]
    missing = [f for f in req_files if not os.path.exists(f)]
    if missing:
        logging.error(f"[{species_name}] Missing files: {missing}")
        return filename
    try:
        sp_fasta_path = os.path.join(LEADING_FASTA_OUT, f"{species_name}_leading_sequences.fasta")
        sp_meta_path = os.path.join(LEADING_FASTA_OUT, f"{species_name}_leading_meta.txt")
        sp_energy_path = os.path.join(LEADING_FASTA_OUT, f"{species_name}_energy.tsv")

        sd_fasta_dict = read_fasta(rRNA_end_sequence_file)
        if not sd_fasta_dict:
            logging.error(f"[{species_name}] SD motif fasta empty")
            return filename
        all_gene_count = count_all_genes(gene_annotation_file)

        leading_gene_count, processed_gene_pairs, leading_regions, leading_rbs_counter = count_leading_genes(
            gene_annotation_file, species_name, GFF3_FOLDER, FASTA_FOLDER, sp_energy_path
        )

        leading_ratio = safe_div(leading_gene_count, all_gene_count)
        strong_sd = leading_rbs_counter.get("Strong SD", 0)
        basic_sd = leading_rbs_counter.get("Basic SD", 0)
        total_sd = strong_sd + basic_sd
        sd_ratio = safe_div(total_sd, leading_gene_count)
        str_sd_ratio = safe_div(strong_sd, leading_gene_count)

        res_dict = {
            'Species': species_name,
            'All Gene Count': all_gene_count,
            'Leading Gene Count': leading_gene_count,
            'Leading Gene Ratio': leading_ratio,
            'Leading Gene with SD': total_sd,
            'Leading Gene with SD%': sd_ratio,
            'Strong SD Count (Leading)': strong_sd,
            'Leading Gene with str_SD%': str_sd_ratio,
            'Basic SD Count (Leading)': basic_sd
        }
        append_single_row_to_csv(LEADING_STAT_CSV, res_dict)

        if leading_regions:
            with open(sp_fasta_path, "w", encoding="utf-8") as fa_out, open(sp_meta_path, "w", encoding="utf-8") as meta_out:
                meta_out.write("gene_id\tfull_dna_seq\tgc_ratio\txi\ttheta\tmfe_energy\trbs_strength\n")
                for gid, data in leading_regions.items():
                    dna_seq = data["RNA_binding_region"][0]
                    fa_out.write(f">{gid}\n{dna_seq}\n")
                    gc_ratio = calculate_gc_content(dna_seq) / 100
                    xi, theta, mfe = data["RNA_binding_region"][2]
                    rbs_strength = data["RNA_binding_region"][3]
                    if mfe is None:
                        meta_out.write(f"{gid}\t{dna_seq}\t{gc_ratio:.4f}\tNA\tNA\tNA\t{rbs_strength}\n")
                    else:
                        meta_out.write(f"{gid}\t{dna_seq}\t{gc_ratio:.4f}\t{xi}\t{theta}\t{mfe}\t{rbs_strength}\n")
            logging.info(f"[{species_name}] Leading fasta + meta + energy saved")
        else:
            logging.warning(f"[{species_name}] No leading genes found")

        logging.info(f"[{species_name}] Finish, Leading count={leading_gene_count}")
        mark_species_done(LEADING_DONE_RECORD, species_name)
        return None
    except Exception as e:
        logging.error(f"[{species_name}] Process failed: {str(e)}", exc_info=True)
        return filename

if __name__ == '__main__':
    done_set = load_finished_species(LEADING_DONE_RECORD)
    all_gff = [f for f in os.listdir(GFF3_FOLDER) if f.endswith(".gff3")]

    todo_gff = []
    for gff_name in all_gff:
        sp_name = gff_name.rsplit(".", 1)[0]
        if sp_name not in done_set:
            todo_gff.append(gff_name)

    print(f"Total gff3 files: {len(all_gff)} | Finished species: {len(done_set)} | Need to run: {len(todo_gff)}")
    if len(todo_gff) == 0:
        print("All Leading gene species already processed, exit directly.")
        exit()

    os.makedirs(LEADING_FASTA_OUT, exist_ok=True)
    unprocessed = set()

    with ThreadPoolExecutor(max_workers=LEADING_MAX_WORKERS) as executor:
        futures = [executor.submit(process_file, fn) for fn in todo_gff]
        for fut in futures:
            err = fut.result()
            if err is not None:
                unprocessed.add(err)

    if unprocessed:
        with open(ERROR_REPORT, "w", encoding="utf-8") as f:
            f.write("Unprocessed gff files:\n")
            for fname in sorted(unprocessed):
                f.write(f"{fname}\n")

    if os.path.exists(LEADING_STAT_CSV):
        df = pd.read_csv(LEADING_STAT_CSV)
        df.to_excel(LEADING_STAT_EXCEL, index=False)
        print(f"Stat excel saved: {LEADING_STAT_EXCEL}")
    else:
        print("No csv stat file found, Excel not generated.")

    print("All Leading gene analysis done.")
    print(f"Fasta & energy output dir: {LEADING_FASTA_OUT}")
