#!/usr/bin/env python3

import click
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped
from Bio import motifs

ALPHABET = Gapped(IUPAC.ambiguous_dna, gap_char='-')
UNAMBIGUOUS_ALPHABET = Gapped(IUPAC.unambiguous_dna, gap_char='-')


def cleanup_alignment(alignment):
    alignment_modified = MultipleSeqAlignment([])

    for rec in alignment:
        if rec.id == '#=GC_RF':
            break
        alignment_modified.append(rec)

    for rec in alignment_modified:
        rec.seq = rec.seq.upper()
        rec.seq = Seq(str(rec.seq).replace('.', '-'), ALPHABET)

    return alignment_modified


def determine_gap_positions(motif, is_unchangeable):
    pwm = motif.counts.normalize()
    gap_counts = pwm['-']
    gap_positions = []

    threshold = 1.0 if is_unchangeable else 0.95
    for i in range(len(gap_counts)):
        if gap_counts[i] >= threshold:
            gap_positions.append(i)

    return gap_positions


def remove_gaps_from_alignment(alignment, species):
    alignment_for_species = MultipleSeqAlignment([rec for rec in alignment if rec.id in species], ALPHABET)

    whole_motif = motifs.create([rec.seq for rec in alignment], ALPHABET)
    species_motif = motifs.create([rec.seq for rec in alignment_for_species], ALPHABET)

    gap_positions_species = determine_gap_positions(species_motif, True)
    gap_positions_whole = determine_gap_positions(whole_motif, False)
    positions_to_remove = set(gap_positions_species).intersection(gap_positions_whole)

    for rec in alignment:
        string_wo_gaps = ''.join([rec.seq[i] for i in range(len(rec.seq)) if i not in positions_to_remove])
        rec.seq = Seq(string_wo_gaps, ALPHABET)

    return alignment


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--rdb-alignment', '-a', type=click.File('r'), help='Alignment in FASTA format')
@click.option('--species', '-sp', type=click.File('r'),
              help='File with taxon representatives each on its own line')
@click.option('--output', '-o', type=click.File('w'))
def foo(rdb_alignment, species, output):
    alignment = AlignIO.read(rdb_alignment, format='fasta')
    alignment_wo_extra_characters = cleanup_alignment(alignment)
    species = [line.strip() for line in species]
    alignment_wo_gaps = remove_gaps_from_alignment(alignment_wo_extra_characters, species)
    AlignIO.write(alignment_wo_gaps, output, format='fasta')

if __name__ == '__main__':
    foo()
