from collections import namedtuple
from alphabet_manipulations import ALPHABET, get_normalized_counts
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio import motifs


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


def remove_gaps(alignment, species):
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


def delete_gaps_from_consensus(alignment):
    alignment_motif = motifs.create([rec.seq for rec in alignment], ALPHABET)
    gap_positions = []
    for i, letter in enumerate(str(alignment_motif.consensus)):
        if letter == '-':
            gap_positions.append(i)

    for rec in alignment:
        seq_string_wo_gaps = ''.join([rec.seq[i] for i in range(len(rec.seq)) if i not in gap_positions])
        rec.seq = Seq(seq_string_wo_gaps, ALPHABET)

    return alignment


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


def create_motifs_from_alignment(alignment, k=20):
    MotifWithInfo = namedtuple('MotifWithPos', ['motif', 'start', 'end', 'ncounts'])
    all_motifs = []
    s = 0
    while s <= (len(alignment[0]) - k):
        motif = motifs.create([rec.seq for rec in alignment[:, s:s + k]], ALPHABET)
        all_motifs.append(MotifWithInfo(motif=motif, start=s, end=s + k, ncounts=get_normalized_counts(motif)))
        s += 1

    return all_motifs
