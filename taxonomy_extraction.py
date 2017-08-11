import os.path as path

__SPECIES2TAXONOMY = None
__DB_PATH = path.join(path.dirname(__file__), 'HITdb_v1.00', 'HITdb_taxonomy_qiime.txt')


def get_species2taxonomy():
    global __SPECIES2TAXONOMY
    if __SPECIES2TAXONOMY is None:
        with open(__DB_PATH) as inf:
            __SPECIES2TAXONOMY = {}
            for line in inf:
                species, taxonomy = line.strip().split(sep="\t")
                __SPECIES2TAXONOMY[species] = taxonomy
    return __SPECIES2TAXONOMY


def get_taxonomy(species):
    return get_species2taxonomy()[species]
