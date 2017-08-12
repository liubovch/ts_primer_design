from itertools import product
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.motifs.matrix import FrequencyPositionMatrix
from Bio import motifs


ALPHABET = Gapped(IUPAC.ambiguous_dna, gap_char='-')
UNAMBIGUOUS_ALPHABET = Gapped(IUPAC.unambiguous_dna, gap_char='-')


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
        if nucl_counts[0] > 0.9:
            key = nucleotides[0]
        elif nucl_counts[1] < 0.1 and sum(nucl_counts[1:]) < 0.25:
            key = nucleotides[0]
        elif nucl_counts[2] < 0.1 and nucl_counts[1] >= 0.1:
            key = ''.join(sorted(nucleotides[:2]))
        elif nucl_counts[3] < 0.1 and nucl_counts[2] >= 0.1:
            key = ''.join(sorted(nucleotides[:3]))
        else:
            key = "ACGT"
        nucleotide = DEGENERATE_NUCL_MAP.get(key, key)
        sequence += nucleotide
    return sequence


def expand_degenerate_sequence(consensus):
    options = [DEGENERATE_NUCL_MAP_INV[x] for x in consensus]
    result = [''.join(x) for x in product(*options)]
    return result


def expand_degenerate_primers(sequences):
    variants = []
    for i_r, rec in enumerate(sequences):
        for i_v, variant in enumerate(expand_degenerate_sequence(rec.seq)):
            variants.append(SeqRecord(Seq(variant), id=f'{rec.id}:{i_r}:{i_v}', description='Primer_to_16S_rRNA'))

    return variants


def get_normalized_counts(motif):
    motif_counts = {}
    for n, counts in motif.counts.items():
        if n in UNAMBIGUOUS_ALPHABET.letters:
            motif_counts[n] = counts
    motif_counts = FrequencyPositionMatrix(UNAMBIGUOUS_ALPHABET, motif_counts)
    normalized_counts = motif_counts.normalize()
    return normalized_counts


def search_positions_w_gaps(primer, sequence):
    primer_id = primer.id.split(':')[0]
    if primer_id in ['R', 'Z_R']:
        i = sequence.ungap().find(primer.seq.reverse_complement())
    else:
        i = sequence.ungap().find(primer.seq)

    if i == -1:
        return None
    else:
        i_start = [i for i, nc in enumerate(sequence) if nc not in ALPHABET.gap_char][i]
        i_end = [i for i, nc in enumerate(sequence) if nc not in ALPHABET.gap_char][i + len(primer.seq) - 1]
        return [i_start, (i_end + 1)]


def create_pwm(alignment, start, end, direction):
    motif = motifs.create([s.seq for s in alignment[:, start:end]], ALPHABET)
    if direction in ['R', 'Z_R']:
        motif = motif.reverse_complement()
    normalized_counts = get_normalized_counts(motif)

    return motif.consensus, normalized_counts
