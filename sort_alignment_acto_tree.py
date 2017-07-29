#!/usr/bin/env python3

import click
from Bio import Phylo
from Bio import AlignIO


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--alignment', '-a', type=click.File('r'), help='Alignment in FASTA format')
@click.option('--species', '-sp', type=click.File('r'),
              help='File with taxon representatives each on its own line')
@click.option('--tree', '-tr', type=click.File('r'),
              help='Phylogenetic tree for species in alignment in NEWICK format')
@click.option('--output', '-o', type=click.File('w'))
def sort_alignment_acto_tree(alignment, species, tree, output):
    tree = Phylo.read(tree, format='newick')
    leaves = tree.get_terminals()

    species = [line.strip() for line in species]
    leaves.sort(key=lambda x: tree.distance(species[0], x))
    ordered_species = species + [item.name for item in leaves if item.name not in species]

    alignment = AlignIO.read(alignment, format='fasta')
    alignment.sort(key=lambda x: ordered_species.index(x.id))
    AlignIO.write(alignment, output, format='fasta')

if __name__ == '__main__':
    sort_alignment_acto_tree()
