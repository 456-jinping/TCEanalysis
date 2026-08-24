# -*- coding: utf-8 -*-
import os
import csv
import logging
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from utils import safe_div, read_fasta, count_all_genes, global_print_lock, load_finished_species, mark_species_done, calculate_gc_content
from core_UTNI import count_UTNI_genes
from config import FASTA_FOLDER, GFF3_FOLDER, RRNA_MOTIF_FOLDER, UTNI_FASTA_OUT, ERROR_REPORT, UTNI_STAT_EXCEL, SD_MOTIF_SUFFIX, UTNI_MAX_WORKERS, UTNI_DONE_RECORD
import threading
csv_write_lock = threading.Lock()
energy_write_lock = threading.Lock()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
results = []

UTNI_STAT_CSV = UTNI_STAT_EXCEL.replace('.xlsx', '.csv')

def append_single_row_to_csv(out_csv, row_dict):
    columns = ['Species', 'All Gene Count', 'UTNI TCEs Count', 'UTNI TCEs Ratio',
               'UTNI TCEs with SD', 'UTNI TCEs with SD%', 'Strong SD Count (UTNI)',
               'UTNI TCEs with str_SD%', 'Basic SD Count (UTNI)']
    with csv_write_lock:
        file_exists = os.path.isfile(out_csv)
        with open(out_csv, 'a', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=columns)
            if not file_exists:
                writer.writeheader()
            writer.writerow(row_dict)
            
def process_file(filename, fasta_folder, gff3_folder, rRNA_folder, utni_fasta_outdir):
    global results
    if not filename.endswith('.gff3'):
        return filename
    species_name = filename.rsplit('.', 1)[0]
    genome_fasta = os.path.join(fasta_folder, f'{species_name}.fasta')
    gene_annotation_file = os.path.join(gff3_folder, f'{species_name}.gff3')
    rRNA_end_sequence_file = os.path.join(rRNA_folder, f"{species_name}{SD_MOTIF_SUFFIX}")
    req_files = [genome_fasta, gene_annotation_file, rRNA_end_sequence_file]
    missing = [f for f in req_files if not os.path.exists(f)]
    if missing:
        logging.error(f"[{species_name}] Missing files: {missing}")
        return filename
    try:
        os.makedirs(utni_fasta_outdir, exist_ok=True)
        sd_fasta_dict = read_fasta(rRNA_end_sequence_file)
        if not sd_fasta_dict:
            logging.error(f"[{species_name}] SD motif fasta empty")
            return filename
        rRNA_end_sequence = list(sd_fasta_dict.values())[0]
        all_gene_count = count_all_genes(gene_annotation_file)

        UTNI_gene_count, processed_gene_pairs, UTNI_genes, UTNI_rbs_strength_counter, seq_dict = count_UTNI_genes(
            gene_annotation_file, species_name, rRNA_folder, gff3_folder, fasta_folder
        )

        UTNI_gene_ratio = safe_div(len(UTNI_genes), all_gene_count)
        UTNI_strong_sd_count = UTNI_rbs_strength_counter.get("Strong SD", 0)
        UTNI_basic_sd_count = UTNI_rbs_strength_counter.get("Basic SD", 0)
        UTNI_Gene_with_SD = UTNI_strong_sd_count + UTNI_basic_sd_count
        UTNI_Gene_with_SD_Ratio = safe_div(UTNI_Gene_with_SD, UTNI_gene_count)
        UTNI_Gene_with_strSD_Ratio = safe_div(UTNI_strong_sd_count, UTNI_gene_count)
        result = {
            'Species': species_name,
            'All Gene Count': all_gene_count,
            'UTNI TCEs Count': UTNI_gene_count,
            'UTNI TCEs Ratio': UTNI_gene_ratio,
            'UTNI TCEs with SD': UTNI_Gene_with_SD,
            'UTNI TCEs with SD%': UTNI_Gene_with_SD_Ratio,
            'Strong SD Count (UTNI)': UTNI_strong_sd_count,
            'UTNI TCEs with str_SD%': UTNI_Gene_with_strSD_Ratio,
            'Basic SD Count (UTNI)': UTNI_basic_sd_count
        }
        with global_print_lock:
            results.append(result)
            append_single_row_to_csv(UTNI_STAT_CSV, result)

        if seq_dict:
            sp_fasta_path = os.path.join(utni_fasta_outdir, f"{species_name}.fasta")
            sp_meta_path = os.path.join(utni_fasta_outdir, f"{species_name}.txt")
            with open(sp_fasta_path, "w", encoding="utf-8") as fa_out, open(sp_meta_path, "w", encoding="utf-8") as meta_out:
                meta_out.write("gene_id\tfull_dna_seq\tgc_ratio\txi\ttheta\tmfe_energy\trbs_strength\n")
                for gid, data in seq_dict.items():
                    dna_seq = data["full_RNA_binding_region"][0]
                    fa_out.write(f">{gid}\n{dna_seq}\n")
                    gc_ratio = calculate_gc_content(dna_seq) / 100
                    xi, theta, mfe = data["full_RNA_binding_region"][2]
                    rbs_strength = data["full_RNA_binding_region"][3]
                    if mfe is None:
                        meta_out.write(f"{gid}\t{dna_seq}\t{gc_ratio:.4f}\tNA\tNA\tNA\t{rbs_strength}\n")
                    else:
                        meta_out.write(f"{gid}\t{dna_seq}\t{gc_ratio:.4f}\t{xi}\t{theta}\t{mfe}\t{rbs_strength}\n")
            logging.info(f"[{species_name}] FASTA + full meta file saved to {utni_fasta_outdir}")
        else:
            logging.warning(f"[{species_name}] No UTNI genes found, skip fasta & meta")

        logging.info(f"[{species_name}] Process complete, UTNI count={UTNI_gene_count}")
        mark_species_done(UTNI_DONE_RECORD, species_name)
        return None
    except Exception as e:
        logging.error(f"[{species_name}] Process failed: {str(e)}", exc_info=True)
        return filename

if __name__ == '__main__':
    results.clear()
    done_set = load_finished_species(UTNI_DONE_RECORD)
    all_species = [f for f in os.listdir(GFF3_FOLDER) if f.endswith(".gff3")]

    todo_gff = []
    for gff_name in all_species:
        sp_name = gff_name.rsplit(".", 1)[0]
        if sp_name not in done_set:
            todo_gff.append(gff_name)

    print(f"Total gff3 files: {len(all_species)} | Finished species: {len(done_set)} | Need to run: {len(todo_gff)}")
    if len(todo_gff) == 0:
        print("All UTNI species already processed, exit directly.")
        exit()

    unprocessed_species = set()
    with ThreadPoolExecutor(max_workers=UTNI_MAX_WORKERS) as executor:
        futures = [
            executor.submit(
                process_file,
                fn,
                FASTA_FOLDER,
                GFF3_FOLDER,
                RRNA_MOTIF_FOLDER,
                UTNI_FASTA_OUT
            )
            for fn in todo_gff
        ]
        for fut in futures:
            res = fut.result()
            if res is not None:
                unprocessed_species.add(res)
    if unprocessed_species:
        with open(ERROR_REPORT, 'w', encoding="utf-8") as f:
            f.write("These species unprocessed / missing file: \n")
            for sp in sorted(unprocessed_species):
                f.write(f"{sp}\n")
    if results:
        df = pd.DataFrame(results)
        df.to_excel(UTNI_STAT_EXCEL, index=False)
        print(f"Final Excel saved to {UTNI_STAT_EXCEL}")
    else:
        print("No results to save, Excel not generated.")
        
#   df = pd.DataFrame(results)
#   df.to_excel(UTNI_STAT_EXCEL, index=False)
    print("All UTNI analysis finished, excel generated.")
    print(f"UTNI gene fasta files saved in dir: {UTNI_FASTA_OUT}")