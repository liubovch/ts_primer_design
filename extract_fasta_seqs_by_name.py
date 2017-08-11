#!/usr/bin/env python3

import click
import os.path as path
from Bio import SeqIO


__SEQUENCES_PATH = path.join(path.dirname(__file__), 'HITdb_v1.00', 'HITdb_sequences.fasta')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--species', '-sp', type=click.File('r'),
              help='File with taxon representatives each on its own line', required=True)
@click.option('--output', '-o', type=click.File('w'), required=True)
def extract_fasta_seqs_by_name(species, output):
    species = [line.strip() for line in species]
    for rec in SeqIO.parse(__SEQUENCES_PATH, 'fasta'):
        if rec.id in species:
            print('>' + rec.id, file=output)
            print(str(rec.seq), file=output)


if __name__ == '__main__':
    extract_fasta_seqs_by_name()
