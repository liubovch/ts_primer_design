#!/usr/bin/env python3

import click
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import motifs
from Bio.motifs.matrix import FrequencyPositionMatrix

from remove_gaps_from_alignment import ALPHABET, UNAMBIGUOUS_ALPHABET


def remove_species_from_alignment(alignment_to_modify, species):
    alignment_modified = MultipleSeqAlignment([], alphabet=ALPHABET)
    for rec in alignment_to_modify:
        if rec.id not in species:
            alignment_modified.append(rec)

    return alignment_modified


def do_alignment_for_taxon(whole_alignment, species):
    taxon_alignment = MultipleSeqAlignment([], alphabet=ALPHABET)
    for rec in whole_alignment:
        if rec.id in species:
            taxon_alignment.append(rec)

    return taxon_alignment


def search_positions_w_gaps(primer, sequence):
    if primer.id in ['R', 'Z_R']:
        i = sequence.ungap().find(primer.seq.reverse_complement())
    else:
        i = sequence.ungap().find(primer.seq)
    i_start = [i for i, nc in enumerate(sequence) if nc not in ALPHABET.gap_char][i]
    i_end = [i for i, nc in enumerate(sequence) if nc not in ALPHABET.gap_char][i + len(primer.seq) - 1]

    return [i_start, (i_end + 1)]


def create_pwm(alignment, start, end, direction):
    motif = motifs.create([s.seq for s in alignment[:, start:end]], ALPHABET)
    if direction in ['R', 'Z_R']:
        motif = motif.reverse_complement()
    normalized_counts = get_normalized_counts(motif)

    return motif.consensus, normalized_counts


def get_normalized_counts(motif):
    motif_counts = {}
    for n, counts in motif.counts.items():
        if n in UNAMBIGUOUS_ALPHABET.letters:
            motif_counts[n] = counts
    motif_counts = FrequencyPositionMatrix(UNAMBIGUOUS_ALPHABET, motif_counts)
    normalized_counts = motif_counts.normalize()
    return normalized_counts


class Report:
    def __init__(self):
        self._lines = []

    def _add_line(self, line=''):
        self._lines.append(line)

    def add_comment(self, name, num_taxon, num_wo_taxon):
        self._add_line('# Normalized position-weight matrices for ' + str(num_taxon) + ' records of ' + name + ' and ' +
                       str(num_wo_taxon) + ' records in all HITdb (without ' + name + ')')
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
@click.option('--alignment', '-a', type=click.File('r'), help='Alignment in FASTA format')
@click.option('--species', '-sp', type=click.File('r'), help='File with taxon representatives each on its own line')
@click.option('--primers', '-p', type=click.File('r'),
              help='List of primers (F — forward, R — reverse, Z_F — forward probe, Z_R — reverse probe)'
                   ' in FASTA format')
@click.option('--output', '-o', type=click.File('w'))
def foo(taxon_name, alignment, species, primers, output):
    whole_alignment = AlignIO.read(alignment, format='fasta', alphabet=ALPHABET)
    species = [line.strip() for line in species]

    alignment_wo_taxon = remove_species_from_alignment(whole_alignment, species)
    taxon_alignment = do_alignment_for_taxon(whole_alignment, species)

    seq, = [r.seq for r in taxon_alignment if r.id == species[0]]

    report = Report()
    report.add_comment(taxon_name, len(taxon_alignment), len(alignment_wo_taxon))
    for rec in SeqIO.parse(primers, 'fasta', alphabet=ALPHABET):
        pos_start, pos_end = search_positions_w_gaps(rec, seq)

        taxon_consensus, taxon_pwm = create_pwm(taxon_alignment, pos_start, pos_end, rec.id)
        wo_taxon_consensus, wo_taxon_pwm = create_pwm(alignment_wo_taxon, pos_start, pos_end, rec.id)

        report.add_primer(rec)
        report.add_consensus(taxon_name, taxon_consensus, taxon_pwm)
        report.add_consensus('w/o ' + taxon_name, wo_taxon_consensus, wo_taxon_pwm)

    output.write(str(report))

if __name__ == '__main__':
    foo()
