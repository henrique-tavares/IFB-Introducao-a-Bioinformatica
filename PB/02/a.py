from os import path
import re
from typing import List, Dict, Pattern


def handle_toxins(file_name: str) -> None:
    complement: Dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}

    sequences_3_5: List[str] = []
    re_toxins: Pattern = re.compile(r">(.*)\n([ATCG\n]*)\n?")

    with open(file_name) as toxins, open("toxins_3-5.fna", "w") as toxins_3_5:
        for match in re_toxins.finditer(toxins.read()):
            sequence: str = match.group(2).replace("\n", "")

            sequence_3_5: str = "".join(complement[n] for n in reversed(sequence))
            sequence_3_5 = "\n".join(
                sequence_3_5[i * 70 : (i + 1) * 70] for i in range((len(sequence_3_5) // 70) + 1)
            ).rstrip("\n")

            sequences_3_5.append(f">{match.group(1)}\n{sequence_3_5}")

        toxins_3_5.write("\n".join(sequences_3_5))


if __name__ == "__main__":
    file_name: str = path.join(path.curdir, "toxinsNCBI.fna")

    try:
        handle_toxins(file_name)
    except IOError:
        print("\nNão foi possíel encontrar o arquivo 'toxinsNCBI.fna' no diretório desse programa\n")
