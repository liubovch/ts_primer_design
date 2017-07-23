#!/usr/bin/env python3

import click
from Bio import SeqIO

@click.command()
@click.option('--species-file', '-sp', type=click.File('r'))
@click.option('--sequences-file', '-seqs', type=click.File('r'))
@click.option('--output-file', '-o', type=click.File('w'))
def extract_fasta_seqs_by_name(species_file, sequences_file, output_file):
    species = [line.strip() for line in species_file]
    for rec in SeqIO.parse(sequences_file, 'fasta'):
        if rec.id in species:
            print('>' + rec.id, file=output_file)
            print(str(rec.seq), file=output_file)


if __name__ == '__main__':
    extract_fasta_seqs_by_name()
