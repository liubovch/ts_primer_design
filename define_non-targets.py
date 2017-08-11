#!/usr/bin/env python3

import click
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from search_for_primer_candidates import expand_degenerate_to_variants, run_blast, collect_hits
from taxonomy_extraction import get_taxonomy


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--primers', '-p', type=click.File('r'),
              help='List of primers (F — forward, R — reverse, Z_F — forward probe, Z_R — reverse probe)'
                   ' in FASTA format presented by ambiguous alphabet')
@click.option('--database-name', '-db', type=str, required=True, help='BLAST database name')
@click.option('--output', '-o', type=click.File('w'))
def define_non_targets(primers, database_name, output):
    sequences = list(SeqIO.parse(primers, format='fasta', alphabet=IUPAC.ambiguous_dna))
    variants = []
    for i_r, rec in enumerate(sequences):
        for i_v, variant in enumerate(expand_degenerate_to_variants(rec.seq)):
            variants.append(SeqRecord(Seq(variant), id=f'{rec.id}:{i_r}:{i_v}', description='Primer_to_16S_rRNA'))

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
