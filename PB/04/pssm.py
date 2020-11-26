import numpy as np
from math import log2
import re
from os import path
import beautifultable as btable
from termcolor import colored

from typing import List, Union, Any, Tuple, Dict
from nptyping import NDArray


def generate_pssm(sequences: List[Tuple[str, int]]) -> NDArray[(4, Any), float]:

    pssm: NDArray[(4, Any), float] = np.ndarray((4, len(sequences[0][0])), dtype=float)

    for position in range(pssm.shape[1]):
        nucleotides: Dict[str, int] = {"A": 0, "T": 0, "G": 0, "C": 0}

        for sequence in sequences:
            nucleotides[sequence[0][position]] += sequence[1]

        for nuc_index, nuc in enumerate(nucleotides.keys()):
            pssm[(nuc_index, position)] = nucleotides[nuc] / sum(nucleotides.values()) or np.nan

    return pssm


def normalize_pssm(pssm: NDArray[(4, Any), float]) -> NDArray[(4, Any), float]:
    normalized_pssm = np.ndarray(pssm.shape, dtype=float)

    for i, row in enumerate(pssm):
        mean = sum(item for item in row if not np.isnan(item)) / len(row)
        normalized_pssm[i] = [item / mean for item in row]

    return normalized_pssm


def convert_pssm_to_log2(normalized_pssm: NDArray[(4, Any), float]) -> NDArray[(4, Any), float]:
    log2_converted_pssm = np.ndarray(normalized_pssm.shape, dtype=float)

    for i, row in enumerate(normalized_pssm):
        log2_converted_pssm[i] = [log2(item) if not np.isnan(item) else np.nan for item in row]

    return log2_converted_pssm


def calculate_score(pssm: NDArray[(4, Any), float], sequence: str) -> float:
    nucleotides: Tuple[str, str, str, str] = ("A", "T", "G", "C")
    score: float = 0

    for j, nucleotide in enumerate(sequence):
        score += pssm[(nucleotides.index(nucleotide), j)] or -1000

    return score


def create_scoring_file(file_name: str, sequences: List[Tuple[str, int]], pssm: NDArray[(Any, Any), float]) -> None:

    with open(file_name, "w") as scoring_file:
        scoring_file.write(f'"Sequência","Score"\n')
        scoring_file.write(
            "\n".join(f'"{sequence[0]}","{calculate_score(pssm, sequence[0]):.3f}"' for sequence in sequences)
        )


def handle_motif(file_name: str) -> List[Tuple[str, int]]:
    sequences: List[Tuple[str, int]] = []

    with open(file_name) as motif:
        re_motif = re.compile(r"(\w+)\s*(\d+)\n?")
        sequences = [(match.group(1), int(match.group(2))) for match in re_motif.finditer(motif.read())]

    if any(len(sequences[0][0]) != len(sequence[0]) for sequence in sequences):
        raise ValueError("As sequências possuem tamanhos diferentes")

    return sequences


if __name__ == "__main__":
    try:
        for file_name in (
            path.join(path.curdir, "C.motif"),
            path.join(path.curdir, "D.motif"),
        ):
            sequences: List[Tuple[str, int]] = handle_motif(file_name)

            print(f"\n{file_name}:\n")

            #####################################################################################################
            #                                   Tabela de frequências simples                                   #
            #####################################################################################################

            pssm = generate_pssm(sequences)

            table_pssm = btable.BeautifulTable()
            table_pssm.set_style(btable.STYLE_BOX)

            for row in pssm:
                table_row: List[str] = [
                    colored(f"{item:.3f}", "green") if not np.isnan(item) else colored("---", "red") for item in row
                ]
                overall_frequency = sum(item for item in row if not np.isnan(item)) / len(row)

                table_pssm.rows.append((*table_row, colored(f"{overall_frequency:.3f}", "cyan")))

            table_pssm.columns.header = [str(i + 1) for i in range(len(sequences[0][0]))] + ["Overall\nfrequency"]
            table_pssm.rows.header = ["A", "T", "G", "C"]

            print("Tabela de frequências simples:", table_pssm, sep="\n", end="\n\n")

            #####################################################################################################
            #                                 Tabela de frequências normalizadas                                #
            #####################################################################################################

            normalized_pssm = normalize_pssm(pssm)

            table_normalized_pssm = btable.BeautifulTable()
            table_normalized_pssm.set_style(btable.STYLE_BOX)

            for row in normalized_pssm:
                table_normalized_pssm.rows.append(
                    (colored(f"{item:.3f}", "green") if not np.isnan(item) else colored("---", "red") for item in row)
                )

            table_normalized_pssm.columns.header = [str(i + 1) for i in range(len(sequences[0][0]))]
            table_normalized_pssm.rows.header = ["A", "T", "G", "C"]

            print("Tabela de frequências normalizadas:", table_normalized_pssm, sep="\n", end="\n\n")

            #####################################################################################################
            #               Tabela de frequências normalizadas para escala logarítimica de base 2               #
            #####################################################################################################

            log2_converted_pssm = convert_pssm_to_log2(normalized_pssm)

            table_log2_converted_pssm = btable.BeautifulTable()
            table_log2_converted_pssm.set_style(btable.STYLE_BOX)

            for row in log2_converted_pssm:
                table_log2_converted_pssm.rows.append(
                    (
                        colored("----", "red")
                        if np.isnan(item)
                        else (colored(f"{item:+.3f}", "green") if item >= 0 else colored(f"{item:+.3f}", "yellow"))
                        for item in row
                    )
                )

            table_log2_converted_pssm.columns.header = [str(i + 1) for i in range(len(sequences[0][0]))]
            table_log2_converted_pssm.rows.header = ["A", "T", "G", "C"]

            print(
                "Tabela de frequências normalizadas para escala logarítimica de base 2:",
                table_log2_converted_pssm,
                sep="\n",
                end="\n\n",
            )

            #####################################################################################################

            create_scoring_file(file_name + ".scores", sequences, log2_converted_pssm)

    except IOError:
        print("Não foi possível encontrar os arquivos 'C.motif' ou 'D.motif' no diretório desse programa", end="\n\n")
