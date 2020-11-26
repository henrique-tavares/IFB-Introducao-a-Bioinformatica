from os import path
import re
from typing import Dict, List


def complement(sequence: str) -> str:
    complement: Dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}

    return "".join((complement[nucleotide] for nucleotide in reversed(sequence)))


def get_orf_by_order(sequence: str, order: int) -> str:

    if order < 0:
        sequence = complement(sequence)
        order *= -1

    orf: str = ""

    for index in range(order - 1, len(sequence), 3):
        codon = sequence[index : index + 3]

        if len(codon) == 3:
            orf += codon

    return orf


def handle_fasta(file_name: str) -> str:
    with open(file_name) as fasta:
        re_fasta = re.compile(r">\w+\n([ATCG]+)\n?")
        sequence: str = re_fasta.match(fasta.read()).group(1)

    return sequence


def generate_multifasta(file_name: str, orfs: Dict[int, str]):
    with open(file_name, "w") as multifasta:
        text: List[str] = []

        for key, value in orfs.items():
            if value:
                text.append(f">orf {key:+d}\n{value}")
            else:
                text.append(f">orf {key:+d}")

        multifasta.write("\n\n".join(text))


if __name__ == "__main__":
    file_name = path.join(path.curdir, "gene.fasta")
    sequence = handle_fasta(file_name)

    orfs: Dict[int, str] = {1: "", 2: "", 3: "", -1: "", -2: "", -3: ""}

    for key in orfs.keys():
        orfs[key] = get_orf_by_order(sequence, key)

    generate_multifasta(path.join(path.curdir, "orfs.fasta"), orfs)
