#!/usr/bin/env python3

import click
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from alphabet_manipulations import expand_degenerate_primers
from search_for_primer_candidates import run_blast, collect_hits
from taxonomy_extraction import get_taxonomy


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--primers', '-p', type=click.File('r'),
              help='List of primers (F — forward, R — reverse, Z_F — forward probe, Z_R — reverse probe)'
                   ' in FASTA format presented by ambiguous alphabet', required=True)
@click.option('--database-name', '-db', type=str, help='BLAST database name', required=True)
@click.option('--output', '-o', type=click.File('w'), required=True)
def define_non_targets(primers, database_name, output):
    sequences = list(SeqIO.parse(primers, format='fasta', alphabet=IUPAC.ambiguous_dna))
    variants = expand_degenerate_primers(sequences)

    blast_results = run_blast(variants, database_name)
    hits = collect_hits(sequences, blast_results)

    output.write(f'Number of non-targets for each primer:\n')
    for rec, non_targets in zip(sequences, hits):
        output.write(f'{rec.id} — {len(non_targets)}\n')

    cross_targets = hits[0].intersection(hits[1].intersection(hits[2]))
    output.write(f'\nNumber of non-targets for primer set: {len(cross_targets)}\n')

    for target in cross_targets:
        output.write(f'{target}\t{get_taxonomy(target)}\n')


if __name__ == '__main__':
    define_non_targets()
