#!/usr/bin/env python3

import click
import sys
from itertools import combinations, groupby
from collections import namedtuple
import tempfile
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from alignment_manipulations import cleanup_alignment, delete_gaps_from_consensus, create_motifs_from_alignment
from alphabet_manipulations import get_degenerate_consensus, expand_degenerate_sequence
from taxonomy_extraction import get_taxonomy


def select_low_degenerate_motifs(all_motifs, num_of_degenerate_nucls=0):
    selected_motifs = []
    for rec in all_motifs:
        conf = 0
        for position in range(len(rec.motif)):
            position_counts = rec.ncounts[:, position]
            if any(v > 0.9 for v in position_counts.values()):
                conf += 1
        if conf >= len(rec.motif) - num_of_degenerate_nucls:
            selected_motifs.append(rec)

    return selected_motifs


def create_list_of_sequences_from_motifs(selected_motifs):
    sequences = []
    for i_r, rec in enumerate(selected_motifs):
        for i_v, variant in enumerate(expand_degenerate_sequence(get_degenerate_consensus(rec.ncounts))):
            sequences.append(
                SeqRecord(Seq(variant), id=f'Motif:{i_r}:{i_v}', description=f'{rec.start}:{rec.end}')
            )

    return sequences


def run_blast(selected_sequences, database_name):
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

    return blast_records


def collect_hits(selected_motifs, blast_records):
    hits = [set() for i in range(len(selected_motifs))]
    for rec in blast_records:
        seqrec_id, _ = rec.query.split(' ')
        _, i, _ = seqrec_id.split(':')
        hits[int(i)].update(alignment.title.split()[0] for alignment in rec.alignments)
    return hits


def choose_motifs_with_min_hits(selected_motifs, blast_records, max_hits=10):
    hits = collect_hits(selected_motifs, blast_records)

    motifs_with_min_hits = []
    for rec, contaminants in zip(selected_motifs, hits):
        if len(contaminants) <= max_hits:
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


MotifPair = namedtuple('MotifPair', ['fst', 'snd', 'cross_hits'])


def choose_motifs_with_min_cross_hits(selected_motifs, blast_records, max_cross_hits=30):
    hits = collect_hits(selected_motifs, blast_records)

    MotifRecWithHits = namedtuple('MotifRecWithHits', ['hits', 'rec'])
    hits_motifs = [MotifRecWithHits(hits=h, rec=r) for h, r in zip(hits, selected_motifs) if len(h) < 500]
    print(f'Number of motifs with less than 500 hits: {len(hits_motifs)}\n')

    motif_pairs = []
    num_cross_hits = []
    for m1, m2 in combinations(hits_motifs, 2):
        if 50 <= abs(m1.rec.start - m2.rec.start) <= 300:
            cross_hits = m1.hits.intersection(m2.hits)
            num_cross_hits.append(len(cross_hits))

            if len(cross_hits) <= max_cross_hits:
                motif_pairs.append(MotifPair(m1.rec, m2.rec, cross_hits))

    print(f'100 minimal numbers of cross-hits: {sorted(num_cross_hits)[:100]}\n')
    return motif_pairs


def simplify_by_first(motif_pairs):
    motif_pairs.sort(key=lambda x: x[0])
    for fst, fst_pairs in groupby(motif_pairs, lambda x: x[0]):
        snds = list(snd for _, snd in fst_pairs)
        edges = {x[:-1]: x[1:] for x in snds}
        incoming = set(edges.values())

        def get_joined(x):
            return x if x not in edges else x[0] + get_joined(edges[x])
        for snd in snds:
            if snd[:-1] not in incoming:
                yield (fst, get_joined(snd[:-1]))


def leave_superstring_first(motif_pairs):
    motif_pairs.sort(key=lambda x: x[0])
    for i in range(len(motif_pairs) - 1, 0, -1):
        if motif_pairs[i][0].startswith(motif_pairs[i - 1][0]):
            motif_pairs[i - 1] = (motif_pairs[i][0], motif_pairs[i - 1][1])


def simplify_pairs(motif_pairs):
    motif_pairs = list(simplify_by_first(motif_pairs))
    motif_pairs = [(y, x) for x, y in motif_pairs]
    leave_superstring_first(motif_pairs)
    motif_pairs = list(simplify_by_first(motif_pairs))
    motif_pairs = [(y, x) for x, y in motif_pairs]
    leave_superstring_first(motif_pairs)
    return motif_pairs


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
@click.option('--max-hits', '-mh', default=10, show_default=True)
@click.option('--max-cross-hits', '-mc', default=30, show_default=True)
def search_for_primer_candidates(alignment, database_name, output, length_of_motif,
                                 max_degenerate, method, max_hits, max_cross_hits):
    print('Job command: ', end='')
    print(*sys.argv)

    alignment = AlignIO.read(alignment, format='fasta')
    taxon_alignment_modified = cleanup_alignment(alignment)
    taxon_alignment_wo_gaps = delete_gaps_from_consensus(taxon_alignment_modified)

    taxon_motifs = create_motifs_from_alignment(taxon_alignment_wo_gaps, length_of_motif)
    print(f'{len(taxon_motifs)} taxon motifs have been created')
    selected_motifs = select_low_degenerate_motifs(taxon_motifs, max_degenerate)
    print(f'{len(selected_motifs)} motifs have been selected')
    selected_sequences = create_list_of_sequences_from_motifs(selected_motifs)
    print(f'{len(selected_sequences)} sequences in all')

    blast_records = run_blast(selected_sequences, database_name)

    if method == 'hits':
        best_motifs = choose_motifs_with_min_hits(selected_motifs, blast_records, max_hits)
        result = concatenate_overlapping_motifs(best_motifs)
        print(f'{len(result)} continuous sequences have been written to the output file')
        for seq in result:
            output.write(seq + '\n')
    else:
        motif_pairs = choose_motifs_with_min_cross_hits(selected_motifs, blast_records, max_cross_hits)
        for motif_pair in motif_pairs:
            print(
                get_degenerate_consensus(motif_pair.fst.ncounts),
                get_degenerate_consensus(motif_pair.snd.ncounts),
                len(motif_pair.cross_hits)
            )

            for hit in motif_pair.cross_hits:
                print(hit, get_taxonomy(hit), sep='\t')
            print()
        print(f'{len(motif_pairs)} pairs of motifs have been selected\n')
        motif_seq_pairs = [(get_degenerate_consensus(x.fst.ncounts), get_degenerate_consensus(x.snd.ncounts)) for x in motif_pairs]
        motif_seq_pairs = simplify_pairs(motif_seq_pairs)
        for fst, snd in motif_seq_pairs:
            output.write(f'{fst}\t{snd}\n')


if __name__ == '__main__':
        search_for_primer_candidates()
