


from csv import DictReader
import sys
from rna_voila.voila_log import voila_log
from tqdm import tqdm

def parse_clin_tsv(tsv_path):
    dataset = {}
    with open(tsv_path, 'r') as f:
        num_lines = sum(1 for _ in f if not _.startswith('#')) - 1

    with open(tsv_path, 'r') as fr:
        reader = DictReader((row for row in fr if not row.startswith('#')), delimiter='\t')
        required_fieldnames = ('gene_id', 'ref_exon_start', 'ref_exon_end', 'event_type', 'event_denovo')
        if not all(x in reader.fieldnames for x in required_fieldnames):
            voila_log().critical(f"CLIN controls tsv requires headers {required_fieldnames} which were not all found.")
            print(reader.fieldnames)
            sys.exit(1)
        voila_log().info("Parsing CLIN controls tsv...")
        for row in tqdm(reader, total=num_lines):
            if row['event_denovo'] == 'True':
                if not row['gene_id'] in dataset:
                    dataset[row['gene_id']] = set()
                lsv_id = f"{row['gene_id']}:{row['event_type']}:{row['ref_exon_start']}-{row['ref_exon_end']}".encode()
                dataset[row['gene_id']].add(lsv_id)

    return dataset