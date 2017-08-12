#!/usr/bin/env python3

import click
from Bio import AlignIO
from Bio import SeqIO
from alphabet_manipulations import ALPHABET, expand_degenerate_primers, search_positions_w_gaps, create_pwm
from alignment_manipulations import remove_species_from_alignment, do_alignment_for_taxon


class Report:
    def __init__(self):
        self._lines = []

    def _add_line(self, line=''):
        self._lines.append(line)

    def add_comment(self, name, num_taxon, num_wo_taxon):
        self._add_line(f'# Normalized position-weight matrices for {str(num_taxon)} records of {name} '
                       f'and {str(num_wo_taxon)} records in all HITdb (without {name}')
        self._add_line()

    def add_primer(self, primer):
        self._add_line('> ' + primer.id)
        self._add_line(str(primer.seq))
        self._add_line()

    def add_consensus(self, name, consensus, counts):
        header = 'Ncl | ' + ''.join('{:^{width}}'.format(cons, width=5) for cons in consensus)
        self._add_line('{:^{width}}'.format('Consensus ' + name + ' (%)', width=len(header)))
        self._add_line('-' * len(header))
        self._add_line(header)
        self._add_line('-' * len(header))
        for nucl, nulc_counts in sorted(counts.items()):
            self._add_line(' {}  | '.format(nucl) + ''.join('{:>3.0f}  '.format(count * 100) for count in nulc_counts))
        self._add_line()

    def __str__(self):
        return '\n'.join(self._lines)


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--taxon-name', '-t', type=str, required=True)
@click.option('--alignment', '-a', type=click.File('r'),
              help='Alignment in FASTA format (after removing gaps)', required=True)
@click.option('--species', '-sp', type=click.File('r'),
              help='File with taxon representatives each on its own line', required=True)
@click.option('--primers', '-p', type=click.File('r'),
              help='List of primers (F — forward, R — reverse, Z_F — forward probe, Z_R — reverse probe)'
                   ' in FASTA format presented by ambiguous alphabet', required=True)
@click.option('--output', '-o', type=click.File('w'), required=True)
def foo(taxon_name, alignment, species, primers, output):
    whole_alignment = AlignIO.read(alignment, format='fasta', alphabet=ALPHABET)
    species = [line.strip() for line in species]

    alignment_wo_taxon = remove_species_from_alignment(whole_alignment, species)
    taxon_alignment = do_alignment_for_taxon(whole_alignment, species)

    seq, = [r.seq for r in taxon_alignment if r.id == species[0]]

    report = Report()
    report.add_comment(taxon_name, len(taxon_alignment), len(alignment_wo_taxon))
    sequences = list(SeqIO.parse(primers, format='fasta', alphabet=ALPHABET))
    variants = expand_degenerate_primers(sequences)
    for variant in variants:
        positions = search_positions_w_gaps(variant, seq)
        if positions is None:
            continue

        taxon_consensus, taxon_pwm = create_pwm(taxon_alignment,
                                                positions[0], positions[1], variant.id.split(':')[0])
        wo_taxon_consensus, wo_taxon_pwm = create_pwm(alignment_wo_taxon,
                                                      positions[0], positions[1], variant.id.split(':')[0])

        report.add_primer(sequences[int(variant.id.split(':')[1])])
        report.add_consensus(taxon_name, taxon_consensus, taxon_pwm)
        report.add_consensus('w/o ' + taxon_name, wo_taxon_consensus, wo_taxon_pwm)

    output.write(str(report))

if __name__ == '__main__':
    foo()
