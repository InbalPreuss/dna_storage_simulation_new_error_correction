from typing import Tuple, Sequence, Generator
import itertools


def dna_sequence_generator(sequence_len=12, symbols=('A', 'C', 'G', 'T')) -> Tuple[str]:
    barcodes = itertools.product(symbols, repeat=sequence_len)
    while True:
        try:
            yield next(barcodes)
        except StopIteration:
            return


def chunker(seq: Sequence, size: int) -> Generator:
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


def decimal_to_bits(decimal_number, amount_bits=None):
    if decimal_number < 0:
        raise ValueError("Only non-negative integers are supported.")

    # Convert to binary and remove the '0b' prefix
    binary_representation = bin(decimal_number)[2:]

    # Pad with leading zeros if amount_bits is specified
    if amount_bits is not None:
        binary_representation = binary_representation.zfill(amount_bits)

    return binary_representation


def iterate_over_bit_chunks(binary_info, chunk_size=3):
    """
    Iterates over a binary string in chunks of a specified size.

    :param binary_info: A string containing the binary information.
    :param chunk_size: The size of each chunk to iterate over (default is 3).
    """
    for i in range(0, len(binary_info), chunk_size):
        chunk = binary_info[i:i + chunk_size]
        # Process the chunk (for now, we'll just print it)
        # print(chunk)


def convert_binary_string_to_tuple(binary_with_info_n_redundancy: str) -> tuple[int, ...]:
    return tuple(int(bit) for bit in binary_with_info_n_redundancy)
