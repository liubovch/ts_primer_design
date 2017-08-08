#!/usr/bin/env python3

import click
from Bio import SeqIO

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--species', '-sp', type=click.File('r'), help='File with taxon representatives each on its own line')
@click.option('--sequences', '-sq', type=click.File('r'),
              help='File in FASTA format from which you want to extract sequences')
@click.option('--output', '-o', type=click.File('w'))
def extract_fasta_seqs_by_name(species, sequences, output):
    species = [line.strip() for line in species]
    for rec in SeqIO.parse(sequences, 'fasta'):
        if rec.id in species:
            print('>' + rec.id, file=output)
            print(str(rec.seq), file=output)


if __name__ == '__main__':
    extract_fasta_seqs_by_name()
