#!/usr/bin/env python3

import click
from itertools import product, combinations
from collections import namedtuple
import tempfile
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import motifs
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from create_PWM_for_primers import get_normalized_counts
from remove_gaps_from_alignment import cleanup_alignment
from remove_gaps_from_alignment import ALPHABET
from genera_extraction import get_genus


DEGENERATE_NUCL_MAP = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'AC': 'M',
    'AG': 'R',
    'AT': 'W',
    'CG': 'S',
    'CT': 'Y',
    'GT': 'K',
    'ACG': 'V',
    'ACT': 'H',
    'AGT': 'D',
    'CGT': 'B',
    'ACGT': 'N',
}

DEGENERATE_NUCL_MAP_INV = {v: k for k, v in DEGENERATE_NUCL_MAP.items()}


# MY OWN FUNCTION for calculating degenerate consensus of motif
def get_degenerate_consensus(counts):
    sequence = ''
    counts = counts
    if '-' in counts:
        del counts['-']
    for i in range(counts.length):
        def get(nucleotide):
            return counts[nucleotide][i]

        nucleotides = sorted(counts, key=get, reverse=True)
        nucl_counts = [counts[c][i] for c in nucleotides]
        if nucl_counts[0] > 0.95:
            key = nucleotides[0]
        elif nucl_counts[1] < 0.05 and sum(nucl_counts[1:]) < 0.08:
            key = nucleotides[0]
        elif nucl_counts[2] < 0.03:
            key = ''.join(sorted(nucleotides[:2]))
        elif nucl_counts[3] < 0.03:
            key = ''.join(sorted(nucleotides[:3]))
        else:
            key = "ACGT"
        nucleotide = DEGENERATE_NUCL_MAP.get(key, key)
        sequence += nucleotide
    return sequence


def expand_degenerate_to_variants(consensus):
    options = [DEGENERATE_NUCL_MAP_INV[x] for x in consensus]
    result = [''.join(x) for x in product(*options)]
    return result


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


def create_motifs_from_alignment(alignment, k=20):
    MotifWithInfo = namedtuple('MotifWithPos', ['motif', 'start', 'end', 'ncounts'])
    all_motifs = []
    s = 0
    while s <= (len(alignment[0]) - k):
        motif = motifs.create([rec.seq for rec in alignment[:, s:s + k]], ALPHABET)
        all_motifs.append(MotifWithInfo(motif=motif, start=s, end=s + k, ncounts=get_normalized_counts(motif)))
        s += 1

    return all_motifs


def select_low_degenerate_motifs(all_motifs, num_of_degenerate_nucls=0):
    selected_motifs = []
    for rec in all_motifs:
        conf = 0
        for position in range(len(rec.motif)):
            position_counts = rec.ncounts[:, position]
            if any(v > 0.95 for v in position_counts.values()):
                conf += 1
        if conf >= len(rec.motif) - num_of_degenerate_nucls:
            selected_motifs.append(rec)

    return selected_motifs


def create_list_of_sequences_from_motifs(selected_motifs):
    sequences = []
    for i_r, rec in enumerate(selected_motifs):
        for i_v, variant in enumerate(expand_degenerate_to_variants(get_degenerate_consensus(rec.ncounts))):
            sequences.append(
                SeqRecord(Seq(variant), id=f'Bif_motif:{i_r}:{i_v}', description=f'{rec.start}:{rec.end}')
            )

    return sequences


def collect_hits(selected_motifs, blast_records):
    hits = [set() for i in range(len(selected_motifs))]
    for rec in blast_records:
        seqrec_id, _ = rec.query.split(' ')
        _, i, _ = seqrec_id.split(':')
        hits[int(i)].update(alignment.title.split()[0] for alignment in rec.alignments)
    return hits


def choose_motifs_with_min_hits(selected_motifs, blast_records, min_hits=10):
    hits = collect_hits(selected_motifs, blast_records)

    motifs_with_min_hits = []
    for rec, contaminants in zip(selected_motifs, hits):
        if len(contaminants) <= min_hits:
            motifs_with_min_hits.append(rec)

    return motifs_with_min_hits


def concatenate_overlapping_motifs(best_motifs):
    sequences = []
    for rec in best_motifs:
        if len(sequences) == 0 or get_degenerate_consensus(rec.ncounts)[:-1] != \
                sequences[-1][len(sequences[-1]) - (len(rec.motif) - 1):]:
            sequences.append(get_degenerate_consensus(rec.ncounts))
        elif get_degenerate_consensus(rec.ncounts)[:-1] == \
                sequences[-1][len(sequences[-1]) - (len(rec.motif) - 1):]:
            sequences[-1] = sequences[-1] + get_degenerate_consensus(rec.ncounts)[-1]

    return sequences


def choose_motifs_with_no_cross_contaminants(selected_motifs, blast_records):
    hits = collect_hits(selected_motifs, blast_records)

    MotifRecWithHits = namedtuple('MotifRecWithHits', ['hits', 'rec'])
    hits_motifs = [MotifRecWithHits(hits=h, rec=r) for h, r in zip(hits, selected_motifs) if len(h) < 500]
    print(f"num motifs less 500: {len(hits_motifs)}")

    results = []
    min_cross_contaminants = None
    contaminants = None
    for m1, m2 in combinations(hits_motifs, 2):
        if 100 <= abs(m1.rec.start - m2.rec.start) <= 300:
            cross_contaminants = m1.hits.intersection(m2.hits)
            if min_cross_contaminants is None or len(cross_contaminants) < min_cross_contaminants:
                min_cross_contaminants = len(cross_contaminants)
                contaminants = cross_contaminants

            if len(cross_contaminants) < 30:
                results.append((
                    get_degenerate_consensus(m1.rec.ncounts),
                    get_degenerate_consensus(m2.rec.ncounts),
                    len(cross_contaminants)
                ))
    print(f'min cross contaminants: {min_cross_contaminants}')
    if contaminants is not None:
        for x in contaminants:
            print(x, get_genus(x))

    return results


def concatenate_cross_contaminants(best_motifs):
    second_column = []


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--alignment', '-a', type=click.File('r'), required=True,
              help='Alignment of taxon species in FASTA format')
@click.option('--database-name', '-db', type=str, required=True, help='BLAST database name')
@click.option('--output', '-o', type=click.File('w'), required=True)
@click.option('--length-of-motif', '-l', default=20, show_default=True)
@click.option('--max-degenerate', '-d', default=0, show_default=True,
              help='Maximum number of degenerate nucleotides allowed in a motif')
@click.option('--method', '-m', type=click.Choice(['hits', 'cross-hits']), required=True)
@click.option('--max-contaminants', '-c', default=10, show_default=True)
def search_for_primer_candidates(alignment, database_name, output, length_of_motif,
                                 max_degenerate, method, max_contaminants):
    alignment = AlignIO.read(alignment, format='fasta')
    taxon_alignment_modified = cleanup_alignment(alignment)
    taxon_alignment_wo_gaps = delete_gaps_from_consensus(taxon_alignment_modified)

    taxon_motifs = create_motifs_from_alignment(taxon_alignment_wo_gaps, length_of_motif)
    print(f'{str(len(taxon_motifs))} taxon motifs have been created')
    selected_motifs = select_low_degenerate_motifs(taxon_motifs, max_degenerate)
    print(f'{str(len(selected_motifs))} motifs have been selected')
    selected_sequences = create_list_of_sequences_from_motifs(selected_motifs)
    print(f'{str(len(selected_sequences))} sequences in all')

    with tempfile.NamedTemporaryFile() as seqs_for_blast, tempfile.NamedTemporaryFile() as blast_results:
        SeqIO.write(selected_sequences, seqs_for_blast.name, format='fasta')
        # run BLAST
        blastn_cline = NcbiblastnCommandline(query=seqs_for_blast.name, db=database_name,
                                             task='blastn-short', perc_identity=100, qcov_hsp_perc=100,
                                             outfmt=5, out=blast_results.name)
        print('Running BLAST!')
        stdout, stderr = blastn_cline()
        print('Extracting BLAST records!')
        results = open(blast_results.name)
        blast_records = list(NCBIXML.parse(results))

    if method == 'hits':
        best_motifs = choose_motifs_with_min_hits(selected_motifs, blast_records, max_contaminants)
        result = concatenate_overlapping_motifs(best_motifs)
        print(f'{str(len(result))} continuous sequences have been written to the output file')
        for seq in result:
            output.write(seq + '\n')
    else:
        best_motifs = choose_motifs_with_no_cross_contaminants(selected_motifs, blast_records)
        for x in best_motifs:
            print(*x)
        # result = concatenate_cross_contaminants(best_motifs)




if __name__ == '__main__':
        search_for_primer_candidates()
