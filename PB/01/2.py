from typing import Dict


def gc_content(sequence: str) -> float:
    if any(nucleotide not in "atgc" for nucleotide in sequence):
        raise AttributeError

    nucleotides: Dict[str, int] = {nucleotide: sequence.count(nucleotide) for nucleotide in "atgc"}

    gc_percentage: float = (nucleotides["g"] + nucleotides["c"]) / sum(nucleotides.values())

    return gc_percentage * 100


if __name__ == "__main__":
    sequence: str = input("\nDigite uma sequência de dna: ").strip().lower()

    try:
        print(f"\n{sequence} -> {gc_content(sequence)}%", end="\n\n")
    except AttributeError:
        print("\nSequencia de dna inválida. Tente novamente.", end="\n\n")

