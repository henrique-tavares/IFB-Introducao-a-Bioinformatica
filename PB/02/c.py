from os import path
import re
from typing import List, Dict, Pattern


def handle_m_rna(file_name: str) -> None:
    protein_translation: Dict[str, str] = {
        "UUU": "Phenylalaine",
        "UUC": "Phenylalaine",
        "UUA": "Leucine",
        "UUG": "Leucine",
        "CUU": "Leucine",
        "CUC": "Leucine",
        "CUA": "Leucine",
        "CUG": "Leucine",
        "AUU": "Isoleucine",
        "AUC": "Isoleucine",
        "AUA": "Isoleucine",
        "AUG": "Methionine",
        "GUU": "Valine",
        "GUC": "Valine",
        "GUA": "Valine",
        "GUG": "Valine",
        "UCU": "Serine",
        "UCC": "Serine",
        "UCA": "Serine",
        "UCG": "Serine",
        "CCU": "Proline",
        "CCC": "Proline",
        "CCA": "Proline",
        "CCG": "Proline",
        "ACU": "Threonine",
        "ACC": "Threonine",
        "ACA": "Threonine",
        "ACG": "Threonine",
        "GCU": "Alanine",
        "GCC": "Alanine",
        "GCA": "Alanine",
        "GCG": "Alanine",
        "UAU": "Tyrosine",
        "UAC": "Tyrosine",
        "UAA": "STOP",
        "UAG": "STOP",
        "CAU": "Histidine",
        "CAC": "Histidine",
        "CAA": "Glutamine",
        "CAG": "Glutamine",
        "AAU": "Asparagine",
        "AAC": "Asparagine",
        "AAA": "Lysine",
        "AAG": "Lysine",
        "GAU": "Aspartic acid",
        "GAC": "Aspartic acid",
        "GAA": "Glutamic acid",
        "GAG": "Glutamic acid",
        "UGU": "Cysteine",
        "UGC": "Cysteine",
        "UGA": "STOP",
        "UGG": "Tryptophan",
        "CGU": "Arginine",
        "CGC": "Arginine",
        "CGA": "Arginine",
        "CGG": "Arginine",
        "AGU": "Serine",
        "AGC": "Serine",
        "AGA": "Serine",
        "AGG": "Serine",
        "GGU": "Glycine",
        "GGC": "Glycine",
        "GGA": "Glycine",
        "GGG": "Glycine",
    }

    protein_sequences: List[str] = []

    re_rna: Pattern = re.compile(r">(.*)\n([UCAG\n]*)\n?")

    with open(file_name) as m_rna, open("toxins.faa", "w") as toxins:
        for match in re_rna.finditer(m_rna.read()):
            proteins: List[str] = []

            rna_sequence: str = match.group(2).replace("\n", "")

            for i in range(len(rna_sequence) // 3):
                protein: str = protein_translation[rna_sequence[i * 3 : (i + 1) * 3]]

                if protein == "STOP":
                    break

                proteins.append(protein)

            protein_list: str = ",\n".join(proteins)
            protein_sequences.append(f">{match.group(1)}\n{protein_list}")

        toxins.write("\n".join(protein_sequences))


if __name__ == "__main__":
    file_name: str = path.join(path.curdir, "mRNA_toxins.fna")

    try:
        handle_m_rna(file_name)
    except IOError:
        print("\nNão foi possíel encontrar o arquivo 'mRNA_toxins.fna' no diretório desse programa\n")
