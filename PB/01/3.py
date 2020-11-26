from typing import Set


def dna_substrings(sequence: str) -> int:
    if any(nucleotide not in "atcg" for nucleotide in sequence):
        raise AttributeError

    substrings: Set[str] = {sequence[i : i + 4] for i in range(0, len(sequence) - 3)}

    return len(substrings)


if __name__ == "__main__":
    sequence: str = input("\nDigite uma sequência de dna: ").strip().lower()

    try:
        print(
            f"\nQuantidades de substrings distintas de tamanho 4 em '{sequence}': {dna_substrings(sequence)}",
            end="\n\n",
        )
    except AttributeError:
        print(f"\nSequencia de dna inválida. Tente novamente.", end="\n\n")
