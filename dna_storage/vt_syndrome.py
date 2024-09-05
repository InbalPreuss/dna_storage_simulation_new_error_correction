from typing import Dict, List, Tuple, Any

import dna_storage.vt_syndrome_utils as vtsu

import copy
import math

import numpy as np


class VTSyndrome:
    def __init__(self, n: int, k: int,
                 bits_per_syndrome: int,
                 payload_len: int,
                 payload_redundancy_len: int,
                 binary_to_k_mer_representation: Dict):
        self.binary_to_k_mer_representation = binary_to_k_mer_representation
        self.n = n
        self.k = k
        self.bits_per_syndrome = bits_per_syndrome

        self.vt_syndrome_utils = vtsu.VTSyndromeUtils(n=n, k=k, bits_per_syndrome=bits_per_syndrome)

        self._input_bits = self.closest_power_of_2()
        self.vt_syn_table = self.vt_syndrome_utils.generate_table()

        self.payload_len = payload_len - payload_redundancy_len

    def codeword_to_syndrome(self, codeword: tuple[int, ...]) -> int:
        syndrome_sum = sum(i * bit for i, bit in enumerate(codeword))
        return np.mod(syndrome_sum, self.n)

    def get_syn_from_input_bits(self, bits_tuple: tuple[int, ...]) -> int:
        return self.vt_syndrome_utils.bits_to_decimal(bits_tuple[:self.bits_per_syndrome])

    def get_index_from_input_bits(self, bits_tuple: tuple[int, ...]) -> int:
        return self.vt_syndrome_utils.bits_to_decimal(bits_tuple[self.bits_per_syndrome:])

    def get_comb_codeword(self, bit_tuple: tuple[int, ...]) -> tuple[int, ...]:
        assert len(bit_tuple) == self._input_bits
        # Step  2 -
        syndrome = self.get_syn_from_input_bits(bits_tuple=bit_tuple)
        idx = self.get_index_from_input_bits(bits_tuple=bit_tuple)
        return self.vt_syn_table[syndrome][idx]

    def closest_power_of_2(self) -> int:
        # Calculate n choose k
        binom = math.comb(self.n, self.k)
        # Find the closest power of 2 less than or equal to binom
        closest_power = np.log2(2 ** (binom.bit_length() - 1))
        return int(closest_power)

    def encode_message(self, message: tuple[int, ...]) -> tuple[int, ...]:
        xs_vector = self.get_comb_codeword(bit_tuple=message)
        return xs_vector

    def encode(self, binary_list) -> tuple[list[int], list[Any]]:
        syndrome_array = []
        xs_vector_array = []
        xs_array = []
        binary_array = [tuple(int(char) for char in string) for string in binary_list]
        for item in binary_array[:self.payload_len]:
            syn, xs_vector = self.get_syn_from_input_bits(
                bits_tuple=item), self.encode_message(message=item)
            syndrome_array.append(syn)
            xs_vector_array.append(xs_vector)
            xs_array.append(self.binary_to_k_mer_representation[item])

        return syndrome_array, xs_array

    def decode(self, codeword: tuple[int, ...], syn_output_from_rs: int):
        current_syn = self.codeword_to_syndrome(codeword=codeword)

        if current_syn != syn_output_from_rs or sum(codeword) < self.k:
            # i = np.mod(np.abs(syn_output_from_rs - current_syn), self.n)
            i = np.mod((int(syn_output_from_rs) - current_syn), self.n) # TODO: check if this is correct, does it needs abs or not?
            codeword_list = list(codeword)
            codeword_list[i] = 1
            codeword = codeword_list
        return tuple(codeword)  # return the corrected codeword


if __name__ == '__main__':
    n = 8
    k = 4
    bits_per_syndrome = 3
    vtsynd = VTSyndrome(n, k, bits_per_syndrome)
    vt_syn_table = vtsynd.vt_syn_table
    vtsynd.vt_syndrome_utils.generate_z_list_from_table()
    # xs_vector is the 8 bits where 4 ones and 4 zeros
    # run the code on 6*6 input to create 6 letters
    bit_tuple = (0, 0, 1, 0, 0, 1)

    # ENCODE
    codeword_encoded = vtsynd.encode_message(message=bit_tuple)
    print(f'encoded codeword: {codeword_encoded}, syn={vtsynd.codeword_to_syndrome(codeword=codeword_encoded)}')

    # NOISE
    noisy_codeword = copy.deepcopy(codeword_encoded)
    noisy_codeword_list = list(noisy_codeword)
    noisy_codeword_list[0] = 0
    noisy_codeword = tuple(noisy_codeword_list)
    print(f'noisy codeword: {noisy_codeword}, syn={vtsynd.codeword_to_syndrome(codeword=noisy_codeword)}')

    # DECODE
    syn_output_from_rs = vtsynd.codeword_to_syndrome(codeword=noisy_codeword)  # not really rs output
    codeword_decoded = vtsynd.decode(codeword=noisy_codeword,
                                     syn_output_from_rs=vtsynd.codeword_to_syndrome(codeword=codeword_encoded))
    print(f'decoded codeword: {codeword_decoded}, syn={vtsynd.codeword_to_syndrome(codeword=codeword_decoded)}')

    # print(syn)
    vtsynd.get_comb_codeword((0, 0, 1, 0, 0, 1))
    vtsynd.get_comb_codeword((0, 0, 1, 0, 1, 1))
    vtsynd.get_comb_codeword((0, 1, 1, 0, 0, 0))
    vtsynd.get_comb_codeword((0, 0, 1, 0, 0, 1))

    x = 4
