from os import path
import re
import numpy as np
import beautifultable as btable
from termcolor import colored

from nptyping import NDArray
from typing import Dict, Any, Union, List, Tuple, Callable


class PairwiseGlobal:
    scoring: Dict[str, int] = {"match": 1, "mismatch": -1, "gap": -2}

    def __init__(self, seq1: str, seq2: str):
        self.__seq1 = seq1
        self.__seq2 = seq2
        self.__matriz = self.calcula_matriz()
        self.__max_score = self.matriz[-1, -1]
        self.__alinhamentos = self.calcula_alinhamentos()

    @property
    def seq1(self) -> str:
        return self.__seq1

    @seq1.setter
    def seq1(self, seq: str):
        self.__seq1 = seq

    @property
    def seq2(self) -> str:
        return self.__seq2

    @seq2.setter
    def seq2(self, seq: str):
        self.__seq2 = seq

    @property
    def matriz(self) -> NDArray[(Any, Any), int]:
        return self.__matriz

    @property
    def max_score(self) -> int:
        return self.__max_score

    @property
    def alinhamentos(self) -> List[List[str]]:
        return self.__alinhamentos

    def match(self, nuc1, nuc2) -> int:
        return self.scoring["match"] if nuc1 == nuc2 else self.scoring["mismatch"]

    def calcula_matriz(self) -> NDArray[(Any, Any), int]:
        matriz: NDArray[(Any, Any), int] = np.zeros((len(self.seq1) + 1, len(self.seq2) + 1), dtype=int)

        for (i, j), _ in np.ndenumerate(matriz):
            if 0 in (i, j):
                matriz[(i, j)] = -2 * max(i, j)
                continue

            matriz[(i, j)] = max(
                matriz[(i - 1, j)] + self.scoring["gap"],
                matriz[(i - 1, j - 1)] + self.match(self.seq1[i - 1], self.seq2[j - 1]),
                matriz[(i, j - 1)] + self.scoring["gap"],
            )

        return matriz

    def calcula_alinhamentos(
        self, bounds: Union[Tuple[int, int], None] = None, alinhamento_base: List[str] = ["", ""],
    ) -> List[List[str]]:
        def safe_scoring(func: Callable[[int, int], int]):
            def wrapper(i: int, j: int) -> Union[int, None]:
                if i < 0 or j < 0:
                    return None
                else:
                    return func(i, j)

            return wrapper

        @safe_scoring
        def up_scoring(i: int, j: int) -> int:
            return self.matriz[(i - 1, j)] + self.scoring["gap"]

        @safe_scoring
        def left_scoring(i: int, j: int) -> int:
            return self.matriz[(i, j - 1)] + self.scoring["gap"]

        @safe_scoring
        def diagonal_scoring(i: int, j: int) -> int:
            return self.matriz[(i - 1, j - 1)] + self.match(self.seq1[i - 1], self.seq2[j - 1])

        alinhamentos: List[List[str]] = []
        alinhamento: List[str] = alinhamento_base

        i: int = bounds[0] if bounds else self.matriz.shape[0] - 1
        j: int = bounds[1] if bounds else self.matriz.shape[1] - 1

        while (i, j) != (0, 0):
            up: Union[int, None] = up_scoring(i, j)
            left: Union[int, None] = left_scoring(i, j)
            diagonal: Union[int, None] = diagonal_scoring(i, j)

            direction = max([score for score in (up, left, diagonal) if score is not None])

            if direction == diagonal:
                if direction == left:
                    alinhamentos += self.calcula_alinhamentos(
                        (i, j - 1), ["-" + alinhamento[0], self.seq2[j - 1] + alinhamento[1]]
                    )

                if direction == up:
                    alinhamentos += self.calcula_alinhamentos(
                        (i - 1, j), [self.seq1[i - 1] + alinhamento[0], "-" + alinhamento[1]]
                    )

                alinhamento[0] = self.seq1[i - 1] + alinhamento[0]
                alinhamento[1] = self.seq2[j - 1] + alinhamento[1]
                i -= 1
                j -= 1

            elif direction == up:
                if direction == left:
                    alinhamentos += self.calcula_alinhamentos(
                        (i, j - 1), ["-" + alinhamento[0], self.seq2[j - 1] + alinhamento[1]]
                    )

                alinhamento[0] = self.seq1[i - 1] + alinhamento[0]
                alinhamento[1] = "-" + alinhamento[1]
                i -= 1

            else:
                alinhamento[0] = "-" + alinhamento[0]
                alinhamento[1] = self.seq2[j - 1] + alinhamento[1]
                j -= 1

        alinhamentos.insert(0, alinhamento)

        return alinhamentos

    @staticmethod
    def handle_fasta(file_name: str) -> List[str]:
        sequences: List[str] = []

        with open(file_name) as fasta:
            re_fasta = re.compile(r">(.*)([ATCG\n]*)\n?")
            sequences = [match.group(2).strip("\n") for match in re_fasta.finditer(fasta.read())]

        return sequences


if __name__ == "__main__":
    file_name: str = path.join(path.curdir, "sequence.fasta")

    try:
        sequences: List[str] = PairwiseGlobal.handle_fasta(file_name)
        seq1: str = sequences[0]
        seq2: str = sequences[1]

        pairwise = PairwiseGlobal(seq1, seq2)

        table_matriz = btable.BeautifulTable()

        for row in pairwise.matriz:
            items: List[str] = []
            for item in row:
                if item >= 0:
                    items.append(colored(str(item), "cyan"))
                else:
                    items.append(colored(str(item), "yellow"))

            table_matriz.rows.append(items)

        table_matriz.columns.header = ["\u03b5"] + [n for n in seq2]
        table_matriz.rows.header = ["\u03b5"] + [n for n in seq1]

        table_matriz.columns.alignment = btable.ALIGN_RIGHT
        table_matriz.set_style(btable.STYLE_BOX_ROUNDED)

        table_alinhamentos = btable.BeautifulTable(maxwidth=100)

        for alinhamento in pairwise.alinhamentos:
            items = ["", ""]
            for i in range(len(alinhamento[0])):
                nuc1 = alinhamento[0][i]
                nuc2 = alinhamento[1][i]

                if pairwise.match(nuc1, nuc2) == pairwise.scoring["match"]:
                    items[0] += colored(nuc1, "green")
                    items[1] += colored(nuc2, "green")
                else:
                    items[0] += colored(nuc1, "red")
                    items[1] += colored(nuc2, "red")

            table_alinhamentos.columns.append(items)

        table_alinhamentos.rows.header = [f"'{seq1}'", f"'{seq2}'"]

        table_alinhamentos.columns.alignment = btable.ALIGN_CENTER
        table_alinhamentos.set_style(btable.STYLE_BOX_ROUNDED)

        print("\nMatriz de Pontuação:")
        print(table_matriz)
        print(f"\nTabela de alinhamentos ótimos, ou seja, com o maior score: {pairwise.max_score}")
        print(table_alinhamentos)

    except IOError:
        print("Não foi possível encontrar o arquivo 'sequence.fasta' no diretório desse programa", end="\n\n")

