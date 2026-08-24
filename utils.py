# -*- coding: utf-8 -*-
import re
import os
import subprocess
import threading
from collections import Counter
from config import EXTEND_UP, EXTEND_DOWN

global_print_lock = threading.Lock()
energy_file_lock = threading.Lock()

# ---------------------- General Tools ----------------------
def safe_div(numerator, denominator):
    if denominator == 0:
        return 0.0
    return float(numerator) / float(denominator)

def read_fasta(fasta_file):
    seq_dict = {}
    with open(fasta_file, 'r') as f:
        seq_id = None
        seq = ''
        for line in f:
            if line.startswith('>'):
                if seq_id:
                    seq_dict[seq_id] = seq
                seq_id = line.strip()[1:]
                seq = ''
            else:
                seq += line.strip()
        if seq_id:
            seq_dict[seq_id] = seq
    return seq_dict

def reverse_complement(seq):
    valid_bases = "ATCG"
    for base in seq:
        if base not in valid_bases:
            raise ValueError(f"Invalid base {base} found in sequence")
    seq = ''.join('A' if base not in valid_bases else base for base in seq)
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'
    }
    return ''.join([complement[base] for base in seq[::-1]])

def calculate_gc_content(seq):
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100

def fetch_transcript_region(raw_chr_seq, start_codon_pos, safe_left, strand, gene_id="", debug=False):
    chr_len = len(raw_chr_seq)
    full_seq = ""
    up_seq = ""
    if strand == '+':
        region_end = start_codon_pos + EXTEND_DOWN - 1
        s = max(0, safe_left)
        e = min(chr_len, region_end)
        full_seq = raw_chr_seq[s:e]
        up_len = start_codon_pos - safe_left - 1
        up_seq = full_seq[:up_len] if up_len > 0 else ""
        if debug:
            print(f"[+] Gene {gene_id}-TCE: seq={full_seq}")
    else:
        s = max(0, start_codon_pos - EXTEND_DOWN)
        e = min(chr_len, safe_left)
        if s >= e:
            return "", ""
        frag = raw_chr_seq[s:e]
        full_seq = reverse_complement(frag)
        up_len = safe_left - start_codon_pos - 1
        up_seq = full_seq[:up_len] if up_len > 0 else ""
        if debug:
                print(f"[-] Gene {gene_id}-TCE:  seq={full_seq}")
    return full_seq, up_seq

def _safe_up_left(start_codon_pos, overlap_gene_end, strand):
    if overlap_gene_end < 0:
        return None
    if strand == '+':
        raw_left = start_codon_pos - EXTEND_UP
        overlap_base = overlap_gene_end - start_codon_pos + 1
        rem = overlap_base % 3
        final_left = raw_left + rem - 1
    else:
        raw_left = start_codon_pos + EXTEND_UP
        overlap_base = raw_left - overlap_gene_end
        rem = overlap_base % 3
        final_left = raw_left - rem
        final_left = min(final_left, overlap_gene_end)
    final_left = max(0, final_left)
    return final_left
    
# ---------------------- RNAhybrid analysis ----------------------
def run_rnahybrid(SD_sequence_file, rRNA_end_sequence_file, xi_values, theta_values):
    try:
        with open(rRNA_end_sequence_file, 'r') as f:
            lines = f.readlines()
            rRNA_end_sequence = ''.join([line.strip() for line in lines[1:]])
        best_energy = None
        best_xi = None
        best_theta = None
        for xi in xi_values:
            for theta in theta_values:
                process = subprocess.Popen(
                    ['RNAhybrid', '-q', SD_sequence_file, '-t', rRNA_end_sequence_file, '-d', f'{xi},{theta}'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
                )
                output, error = process.communicate()
                process.stdout.close()
                process.stderr.close()
                process.wait()
                if error:
                    with global_print_lock:
                        print(f"Error running RNAhybrid with xi={xi}, theta={theta}: {error}")
                    continue
                pattern = re.compile(r'mfe:\s*([+-]?\d+\.?\d*)\s*kcal/mol', re.MULTILINE)
                match = pattern.search(output)
                if match:
                    energy = float(match.group(1))
                    if best_energy is None or energy < best_energy:
                        best_energy = energy
                        best_xi = xi
                        best_theta = theta
#                       with global_print_lock:
#                           print(f"【New best mfe】xi={xi}, theta={theta}, energy={energy}")
                else:
                    with global_print_lock:
                        print(f"No match found in RNAhybrid output with xi={xi}, theta={theta}")
        return best_xi, best_theta, best_energy
    except Exception as e:
        with global_print_lock:
            print(f"RNAhybrid run error: {e}")
        return None, None, None

# ---------------------- GFF analysis tool ----------------------
def read_strand_info(gene_annotation_file):
    strand_dict = {}
    gene_id_pattern = re.compile(r'ID=([^;]+)')
    gene_name_pattern = re.compile(r'Name=([^;]+)')
    with open(gene_annotation_file, 'r') as ann_file:
        for line in ann_file:
            cols = line.strip().split('\t')
            if len(cols) >= 9:
                attributes = cols[8].split(';')
                gene_id = None
                for attr in attributes:
                    if attr.startswith('ID='):
                        gene_id = attr.split('=')[1].strip('"')
                strand = cols[6]
                if gene_id:
                    strand_dict[gene_id] = strand
    return strand_dict

def count_all_genes(gene_annotation_file):
    all_genes = set()
    gene_id_pattern = re.compile(r'ID=([^;]+)')
    gene_name_pattern = re.compile(r'Name=([^;]+)')
    import logging
    try:
        with open(gene_annotation_file, 'r') as ann_file:
            for line in ann_file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                cols = line.split('\t')
                if len(cols) < 9:
                    continue
                feature_type = cols[2]
                if feature_type != "gene":
                    continue
                attr_str = cols[8]
                gene_id = None
                gene_name = None
                id_match = gene_id_pattern.search(attr_str)
                name_match = gene_name_pattern.search(attr_str)
                if id_match:
                    gene_id = id_match.group(1).strip('"')
                if name_match:
                    gene_name = name_match.group(1).strip('"')
                if gene_id:
                    all_genes.add(gene_id)
                elif gene_name:
                    all_genes.add(gene_name)
    except Exception as e:
        logging.error(f"Error counting genes: {str(e)}")
    return len(all_genes)

# ---------------------- UTNI-RNAfold tool ----------------------
def predict_structure_and_energy(seq):
    try:
        process = subprocess.Popen(['RNAfold', '--noPS'], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True)
        output, error = process.communicate(input=seq)
        process.stdout.close()
        process.stderr.close()
        process.wait()
        if error:
            logging.error(f"Error running RNAfold: {error}")
            return None, None
        lines = output.strip().split('\n')
        if len(lines) < 2:
            logging.error("Unexpected output from RNAfold")
            return None, None
        structure = lines[1].split()[0]
        gibbs_free_energy = lines[1].split()[1].strip('()')
        try:
            gibbs_free_energy = float(gibbs_free_energy)
        except ValueError:
            gibbs_free_energy = float('nan')
        return structure, gibbs_free_energy
    except FileNotFoundError:
        logging.error("RNAfold is not installed or not in the system path.")
        return None, None
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        return None, None

def Hairpin_structure(rna_structure):
    for i, char in enumerate(rna_structure):
        if char == '(':
            try:
                close_index = rna_structure.index(')', i + 1)
                return True
            except ValueError:
                return False
    return False

def parse_gene_info(gene_annotation_file):
    gene_info_dict = {}
    start_codon_pos_dict = {}
    try:
        with open(gene_annotation_file, 'r') as ann_file:
            for line in ann_file:
                if line.startswith('#'):
                    continue
                cols = line.strip().split('\t')
                if len(cols) >= 9:
                    feature_type = cols[2]
                    if feature_type != "gene":
                        continue
                    attributes = cols[8].split(';')
                    gene_id = None
                    chromosome = cols[0]
                    start = int(cols[3])
                    end = int(cols[4])
                    strand = cols[6]
                    for attr in attributes:
                        if attr.startswith('ID='):
                            gene_id = attr.split('=')[1].strip('"')
                    if gene_id:
                        if strand == '-':
                            start_codon_pos = end
                        elif strand == '+':
                            start_codon_pos = start
                        gene_info_dict[gene_id] = (start, end, strand, chromosome)
                        start_codon_pos_dict[gene_id] = start_codon_pos
        return gene_info_dict, start_codon_pos_dict
    except FileNotFoundError:
        import logging
        logging.error(f"File not found: {gene_annotation_file}")
        return {}, {}
    except Exception as e:
        import logging
        logging.error(f"Error parsing gene info: {e}")
        return {}, {}
        
def load_finished_species(record_path):
    done_sp = set()
    if os.path.exists(record_path):
        with open(record_path, "r", encoding="utf-8") as f:
            for line in f:
                sp = line.strip()
                if sp:
                    done_sp.add(sp)
    return done_sp

def mark_species_done(record_path, species_name):
    with energy_file_lock:
        with open(record_path, "a", encoding="utf-8") as f:
            f.write(f"{species_name}\n")
            
# ---------------------- safe_up_left analysis ----------------------
def _safe_up_left(start_codon_pos, overlap_gene_end, strand):
    if overlap_gene_end < 0:
        return None
    if strand == '+':
        raw_left = start_codon_pos - EXTEND_UP
        overlap_base = overlap_gene_end - start_codon_pos + 1
        rem = overlap_base % 3
        final_left = raw_left + rem - 1
    else:
        raw_left = start_codon_pos + EXTEND_UP
        overlap_base = raw_left - overlap_gene_end
        rem = overlap_base % 3
        final_left = raw_left - rem
        final_left = min(final_left, overlap_gene_end)
    final_left = max(0, final_left)
    return final_left
