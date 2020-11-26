def fatorial(n: int) -> int:
    if n < 0:
        raise AttributeError
    elif n == 0 or n == 1:
        return 1
    else:
        return n * fatorial(n - 1)


if __name__ == "__main__":
    try:
        n: int = int(input("\nDigite um número para calcular o seu fatorial: "))
        print(f"\n{n}! = {fatorial(n)}", end="\n\n")
    except ValueError:
        print("\nPor favor digite um número valido. Tente novamente.", end="\n\n")
    except AttributeError:
        print("\nNão é possível calcular o fatorial de um número negativo.", end="\n\n")
