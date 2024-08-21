import copy
import math
from itertools import combinations
from typing import List, Tuple

import numpy as np


class VTSyndrome:
    def __init__(self, n: int, k: int):
        self.n = n
        self.k = k

        self._input_bits = self.closest_power_of_2()
        self._vt_syn_table = self.generate_table()

    def create_word(self, ones_positions) -> List[int]:
        bits = [0] * self.n
        for pos in ones_positions:
            bits[pos] = 1
        return bits

    def generate_table(self) -> List[List[List[int]]]:
        table = [[], [], [], [], [], [], [], []]
        for ones_positions in combinations(range(self.n), self.k):
            bits = self.create_word(ones_positions)
            syndrome_index = self.codeword_to_syndrome(codeword=bits)
            if len(table[syndrome_index]) < self.n:
                table[syndrome_index].append(bits)
        return table

    def print_table(self, vt_syn_table: List[List[List[int]]]) -> None:
        for i in range(8):
            print(vt_syn_table[i])
            print(len(vt_syn_table[i]))

    def codeword_to_syndrome(self, codeword: List[int]) -> int:
        syndrome_sum = sum(i * bit for i, bit in enumerate(codeword))
        return np.mod(syndrome_sum, self.n)
    def bits_to_decimal(self, bits: List[int]) -> int:
        syn = 0
        for i, bit in enumerate(reversed(bits)):  # reverse to start from LSB
            syn += bit * (2 ** i)
        return syn

    def get_syn_from_input_bits(self, bits_array: List[int]) -> int:
        split_bit = int(self._input_bits / 2)
        return self.bits_to_decimal(bits_array[:split_bit])

    def get_index_from_input_bits(self, bits_array: List[int]) -> int:
        split_bit = int(self._input_bits / 2)
        return self.bits_to_decimal(bits_array[split_bit:])

    def get_comb_codeword(self, bit_array: List[int]) -> List[int]:
        assert len(bit_array) == self._input_bits
        # Step  2 -
        syndrome = self.get_syn_from_input_bits(bits_array=bit_array)
        idx = self.get_index_from_input_bits(bits_array=bit_array)
        return self._vt_syn_table[syndrome][idx]

    def closest_power_of_2(self) -> int:
        # Calculate n choose k
        binom = math.comb(self.n, self.k)
        # Find the closest power of 2 less than or equal to binom
        closest_power = np.log2(2 ** (binom.bit_length() - 1))
        return int(closest_power)

    def encode(self, message: List[int]) -> list[int]:
        xs_vector = self.get_comb_codeword(bit_array=message)
        return xs_vector

    def decode(self, codeword: List[int], syn_output_from_rs: int):
        current_syn = self.codeword_to_syndrome(codeword=codeword)

        if current_syn != syn_output_from_rs or sum(codeword) < self.k:
            i = np.mod(np.abs(syn_output_from_rs - current_syn), self.n)
            codeword[i] = 1
        return codeword  # return the corrected codeword


if __name__ == '__main__':
    n = 8
    k = 4
    vtsynd = VTSyndrome(n, k)
    vt_syn_table = vtsynd.generate_table()
    # xs_vector is the 8 bits where 4 ones and 4 zeros
    # run the code on 6*6 input to create 6 letters
    bit_array = [0, 0, 1, 0, 0, 1]
    # syn, ind, xs_vector = vtsynd.get_syn_from_input_bits(
    #     bits_array=bit_array), vtsynd.get_index_from_input_bits(
    #     bits_array=bit_array), vtsynd.get_comb_codeword(bit_array=bit_array)
    codeword_encoded = vtsynd.encode(message=bit_array)
    print(f'encoded codeword: {codeword_encoded}, syn={vtsynd.codeword_to_syndrome(codeword=codeword_encoded)}')
    noisy_codeword = copy.deepcopy(codeword_encoded)
    noisy_codeword[0] = 0
    print(f'noisy codeword: {noisy_codeword}, syn={vtsynd.codeword_to_syndrome(codeword=noisy_codeword)}')
    syn_output_from_rs = vtsynd.codeword_to_syndrome(codeword=noisy_codeword) #not really rs output
    codeword_decoded = vtsynd.decode(codeword=noisy_codeword, syn_output_from_rs=vtsynd.codeword_to_syndrome(codeword=codeword_encoded))
    print(f'decoded codeword: {codeword_decoded}, syn={vtsynd.codeword_to_syndrome(codeword=codeword_decoded)}')




    # print(syn)
    vtsynd.get_comb_codeword(vt_syn_table, [0, 0, 1, 0, 0, 1])
    vtsynd.get_comb_codeword(vt_syn_table, [0, 0, 1, 0, 1, 1])
    vtsynd.get_comb_codeword(vt_syn_table, [0, 1, 1, 0, 0, 0])
    vtsynd.get_comb_codeword(vt_syn_table, [0, 0, 1, 0, 0, 1])

    x = 4
