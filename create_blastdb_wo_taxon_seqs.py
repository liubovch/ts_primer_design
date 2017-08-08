#!/usr/bin/env python3

import click
from Bio import SeqIO
import tempfile
import subprocess


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--species', '-sp', type=click.File('r'), help='File with taxon representatives each on its own line')
@click.option('--sequences', '-sq', type=click.File('r'), help='Sequences in FASTA format')
@click.option('--database-name', '-db', help='Name of BLAST database to be created')
def create_blastdb(species, sequences, database_name):
    species = [line.strip() for line in species]

    db_wo_taxon_seqs = []
    for rec in SeqIO.parse(sequences, format='fasta'):
        if rec.id not in species:
            db_wo_taxon_seqs.append(rec)
    with tempfile.NamedTemporaryFile() as tmp:
        SeqIO.write(db_wo_taxon_seqs, tmp.name, format='fasta')
        subprocess.run(['makeblastdb', '-in', tmp.name, '-input_type', 'fasta', '-parse_seqids',
                        '-dbtype', 'nucl', '-out', database_name])

if __name__ == '__main__':
    create_blastdb()
