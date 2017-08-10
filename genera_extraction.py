import os.path as path

__SPECIES2GENUS = None
__DB_PATH = path.join(path.dirname(__file__), 'HITdb_v1.00', 'HITdb_taxonomy_qiime.txt')


def get_species2genus():
    global __SPECIES2GENUS
    if __SPECIES2GENUS is None:
        with open(__DB_PATH) as inf:
            __SPECIES2GENUS = {}
            for line in inf:
                species, taxonomy = line.strip().split(sep="\t")
                *_, genus, rep = taxonomy.split(sep=";")
                if genus[0] == '[':
                    genus = genus[1:-1]
                __SPECIES2GENUS[species] = genus
    return __SPECIES2GENUS


def get_genus(species):
    return get_species2genus()[species]
