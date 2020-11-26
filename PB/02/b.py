from os import path
import re
from typing import List, Dict, Pattern


def handle_toxins_3_5(file_name: str) -> None:
    rna_translation: Dict[str, str] = {"A": "U", "T": "A", "C": "G", "G": "C"}

    re_toxins: Pattern = re.compile(r">(.*)\n([ATCG\n]*)\n?")

    rna_sequences: List[str] = []

    with open(file_name) as toxins_3_5, open("mRNA_toxins.fna", "w") as rna_toxins:
        for match in re_toxins.finditer(toxins_3_5.read()):
            sequence_3_5: str = match.group(2).replace("\n", "")

            m_rna: str = "".join(rna_translation[n] for n in reversed(sequence_3_5))
            m_rna = "\n".join(m_rna[i * 70 : (i + 1) * 70] for i in range((len(m_rna) // 70) + 1)).rstrip("\n")

            rna_sequences.append(f">{match.group(1)}\n{m_rna}")

        rna_toxins.write("\n".join(rna_sequences))


if __name__ == "__main__":
    file_name: str = path.join(path.curdir, "toxins_3-5.fna")

    try:
        handle_toxins_3_5(file_name)
    except IOError:
        print("\nNão foi possíel encontrar o arquivo 'toxins_3-5.fna' no diretório desse programa\n")
