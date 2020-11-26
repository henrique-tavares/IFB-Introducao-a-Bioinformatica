from os import path
from collections import defaultdict
from typing import List, Tuple, DefaultDict
from operator import itemgetter


def handle_fasta_sequence_substring(raw_sequence: str) -> Tuple[str, str, int]:
    sequence_items: List[str] = [item for item in raw_sequence.rstrip("\n").split("\n")]

    header: str = sequence_items[0]
    sequence: str = sequence_items[1]
    substrings: DefaultDict[str, int] = defaultdict(int)

    for i in range(0, len(sequence) - 2):
        substring = sequence[i : i + 3]

        substrings[substring] += 1

    most_frequent_substring = max(substrings.items(), key=itemgetter(1))

    return (header, most_frequent_substring[0], most_frequent_substring[1])


def handle_fasta(file_name: str) -> None:

    with open(file_name) as fasta, open("results.txt", "w") as results:
        raw_sequences: List[str] = fasta.read().lstrip(">").split(">")
        sequences: List[Tuple[str, str, int]] = list(map(handle_fasta_sequence_substring, raw_sequences))

        formatted_sequences: List[str] = list(
            map(lambda sequence: f">{sequence[0]} {sequence[1]} {sequence[2]}", sequences)
        )

        results.write("\n".join(formatted_sequences))


if __name__ == "__main__":
    file_name: str = path.join(path.curdir, "sequence.fasta")

    try:
        handle_fasta(file_name)
    except IOError:
        print("\nNão foi possível encontrar o arquivo 'sequence.fasta' no diretório desse programa", end="\n\n")
