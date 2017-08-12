#!/usr/bin/env python3

import click
from Bio import AlignIO
from alignment_manipulations import cleanup_alignment, remove_gaps


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--alignment', '-a', type=click.File('r'),
              help='Alignment of all sequences from HITdb in FASTA format', required=True)
@click.option('--species', '-sp', type=click.File('r'),
              help='File with taxon representatives each on its own line', required=True)
@click.option('--output', '-o', type=click.File('w'), required=True)
def remove_gaps_from_alignment(alignment, species, output):
    alignment = AlignIO.read(alignment, format='fasta')
    alignment_wo_extra_characters = cleanup_alignment(alignment)
    species = [line.strip() for line in species]
    alignment_wo_gaps = remove_gaps(alignment_wo_extra_characters, species)
    AlignIO.write(alignment_wo_gaps, output, format='fasta')

if __name__ == '__main__':
    remove_gaps_from_alignment()
