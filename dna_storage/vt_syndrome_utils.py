from itertools import combinations
from typing import List, Tuple

import numpy as np

import dna_storage.utils as uts


def convert_to_x_tuple(x_list) -> tuple[str, ...]:
    return tuple(f'X{i + 1}' for i, bit in enumerate(x_list) if bit == 1)


class VTSyndromeUtils:
    def __init__(self, n: int, k: int, bits_per_syndrome: int):
        self.n = n
        self.k = k
        self.bits_per_syndrome = bits_per_syndrome

    def print_table(self, vt_syn_table: List[List[List[int]]]) -> None:
        for i in range(8):
            print(vt_syn_table[i])
            print(len(vt_syn_table[i]))

    def bits_to_decimal(self, bits: tuple[int, ...]) -> int:
        syn = 0
        for i, bit in enumerate(reversed(bits)):  # reverse to start from LSB
            syn += bit * (2 ** i)
        return syn

    def create_word(self, ones_positions) -> tuple[int, ...]:
        bits = [0] * self.n
        for pos in ones_positions:
            bits[pos] = 1
        return tuple(bits)

    def codeword_to_syndrome(self, codeword: tuple[int, ...]) -> int:
        syndrome_sum = sum(i * bit for i, bit in enumerate(codeword))
        return np.mod(syndrome_sum, self.n)

    def generate_table(self) -> List[List[tuple[int, ...]]]:
        table = [[], [], [], [], [], [], [], []]
        for ones_positions in combinations(range(self.n), self.k):
            bits = self.create_word(ones_positions)
            syndrome_index = self.codeword_to_syndrome(codeword=bits)
            if len(table[syndrome_index]) < self.n:
                table[syndrome_index].append(bits)
        return table

    def generate_z_list_from_table(self) -> Tuple[dict, dict, dict, dict, dict, dict, dict]:
        table = self.generate_table()

        z_index = 1
        z_to_k_mer_representative = {}
        k_mer_representative_to_z = {}
        z_to_binary = {}
        binary_to_z = {}
        binary_to_k_mer_representation = {}
        k_mer_representation_to_kmer_vector_representation = {}
        kmer_vector_representation_to_mer_representation = {}
        for syn, xs_list in enumerate(table):
            for index, x_list in enumerate(xs_list):
                z_name = 'Z' + str(z_index)
                x_tuple = convert_to_x_tuple(x_list=x_list)
                z_to_k_mer_representative[z_name] = x_tuple
                k_mer_representative_to_z[x_tuple] = z_name
                kmer_vector_representation_to_mer_representation[x_list] = x_tuple
                k_mer_representation_to_kmer_vector_representation[x_tuple] = x_list

                syn_as_bits = uts.decimal_to_bits(decimal_number=syn, amount_bits=self.bits_per_syndrome)
                index_as_bits = uts.decimal_to_bits(decimal_number=index, amount_bits=self.bits_per_syndrome)
                binary = syn_as_bits + index_as_bits
                binary_tuple = tuple(int(bit) for bit in binary)
                z_to_binary[z_name] = binary_tuple
                binary_to_z[binary_tuple] = z_name
                binary_to_k_mer_representation[binary_tuple] = x_tuple

                z_index += 1

        return (z_to_k_mer_representative,
                k_mer_representative_to_z,
                z_to_binary, binary_to_z,
                binary_to_k_mer_representation,
                k_mer_representation_to_kmer_vector_representation,
                kmer_vector_representation_to_mer_representation)
