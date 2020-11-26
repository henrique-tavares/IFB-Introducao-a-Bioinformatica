from os import path
from typing import List, Tuple


def handle_fasta_tuple(raw_sequence: str) -> Tuple[str, ...]:
    return tuple(item for item in raw_sequence.rstrip("\n").split("\n"))


def handle_fasta(file_name: str) -> List[Tuple[str, ...]]:
    sequences: List[Tuple[str, ...]] = []

    with open(file_name) as fasta:
        raw_sequences: List[str] = fasta.read().lstrip(">").split(">")
        sequences = list(map(handle_fasta_tuple, raw_sequences))

    return sequences


if __name__ == "__main__":
    file_name: str = path.join(path.curdir, "sequence.fasta")

    try:
        for fasta_tuple in handle_fasta(file_name):
            print(fasta_tuple)
    except IOError:
        print("\nNão foi possível encontrar o arquivo 'sequence.fasta' no diretório desse programa", end="\n\n")
