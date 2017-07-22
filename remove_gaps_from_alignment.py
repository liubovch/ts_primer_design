#!/usr/bin/env python3

import click
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped
from Bio import motifs

ALPHABET = Gapped(IUPAC.ambiguous_dna, gap_char='-')


def cleanup_alignment(alignment):
    alignment_modified = MultipleSeqAlignment([])

    for record in alignment:
        if record.id == '#=GC_RF':
            break
        alignment_modified.append(record)

    for record in alignment_modified:
        record.seq = record.seq.upper()
        record.seq = Seq(str(record.seq).replace('.', '-'), ALPHABET)

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
    alignment_for_species = MultipleSeqAlignment([r for r in alignment if r.id in species], ALPHABET)

    whole_motif = motifs.create([r.seq for r in alignment], ALPHABET)
    species_motif = motifs.create([r.seq for r in alignment_for_species], ALPHABET)

    gap_positions_species = determine_gap_positions(species_motif, True)
    gap_positions_whole = determine_gap_positions(whole_motif, False)
    positions_to_remove = set(gap_positions_species).intersection(gap_positions_whole)

    for r in alignment:
        seq_string_wo_gaps = ''.join([r.seq[i] for i in range(len(r.seq)) if i not in positions_to_remove])
        r.seq = Seq(seq_string_wo_gaps, ALPHABET)

    return alignment


@click.command()
@click.option('--rdb-alignment-file', type=click.File('r'))
@click.option('--species-file', type=click.File('r'))
@click.option('--output-file', type=click.File('w'))
def foo(rdb_alignment_file, species_file, output_file):
    alignment = AlignIO.read(rdb_alignment_file, format='fasta')
    alignment_wo_extra_characters = cleanup_alignment(alignment)
    species = [line.strip() for line in species_file]
    alignment_wo_gaps = remove_gaps_from_alignment(alignment_wo_extra_characters, species)
    print(alignment_wo_gaps)
    AlignIO.write(alignment_wo_gaps, output_file, format='fasta')

if __name__ == '__main__':
    foo()
