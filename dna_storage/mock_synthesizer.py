import itertools
from pathlib import Path
import random
from typing import Union, Dict, List

import numpy as np

from dna_storage.utils import chunker


class Synthesizer:
    def __init__(self, input_file: Union[Path, str],
                 results_file: Union[Path, str],
                 synthesis_config: Dict,
                 barcode_total_len: int,
                 subset_size: int,
                 k_mer_representative_to_z: Dict,
                 k_mer_to_dna: Dict,
                 k_mer: int,
                 mode: str):
        self.input_file = input_file
        self.results_file = results_file
        open(self.results_file, 'w').close()
        self.synthesis_config = synthesis_config
        self.barcode_total_len = barcode_total_len
        self.subset_size = subset_size
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.k_mer_to_dna = k_mer_to_dna
        self.k_mer = k_mer
        self.mode = mode

    def synthesize(self):
        if self.mode == 'test':
            np.random.seed(self.synthesis_config['seed'])
            random.seed(self.synthesis_config['seed'])
        with open(self.input_file, 'r', encoding='utf-8') as input_file, open(self.results_file, 'w+', encoding='utf-8') as results_file:
            for line in input_file:
                line_list = line.strip('\n').split(',')
                barcode, payload = line_list[0], line_list[1:]
                x_list = self.get_x_list(payload=payload)
                number_of_nuc = max(1, int(round(np.random.normal(self.synthesis_config['number_of_oligos_per_barcode'],scale=10))))
                x_mat = np.empty([number_of_nuc, self.barcode_total_len], dtype=np.dtype(('U', 5)))
                x_mat[:] = np.array(tuple(barcode))
                for idx, x_tuple in enumerate(x_list, 1):
                    vec = np.random.choice(x_tuple, size=(number_of_nuc,))
                    col = np.array([tuple(self.k_mer_to_dna[k_mer]) for k_mer in vec])
                    x_mat = np.hstack((x_mat, col))

                barcode_list = [''.join(row) for row in x_mat[:, :self.barcode_total_len]]
                payload_list = [''.join(row) for row in x_mat[:, self.barcode_total_len:]]
                barcode_list = self.insertion_deletion_substitution(barcode_list, group_size=1)
                payload_list = self.insertion_deletion_substitution(payload_list, group_size=self.k_mer)
                dna_list = [b + p for b, p in zip(barcode_list, payload_list)]
                results_file.write('\n'.join(dna_list) + '\n')

    def insertion_deletion_substitution(self, dna_list: List[str], group_size: int = 1):
        choose_from = 'ACGT' if group_size == 1 else list(self.k_mer_to_dna.values())
        choose_from_set = set(choose_from)
        bitmap_length = int(len(dna_list[0]) / group_size)
        for row_idx, oligo in enumerate(dna_list):
            deletion = np.random.binomial(1, self.synthesis_config['letter_deletion_error_ratio'], bitmap_length)
            oligo = ''.join([group if deletion[idx] == 0 else '' for idx, group in enumerate(chunker(oligo, group_size))])
            insertion_idx = np.random.binomial(1, self.synthesis_config['letter_insertion_error_ratio'], bitmap_length)
            insertion = [random.choice(choose_from) if i == 1 else '' for i in insertion_idx]
            oligo = ''.join(''.join(x) for x in itertools.zip_longest(chunker(oligo, group_size), insertion, fillvalue=''))
            substitution_idx = np.random.binomial(1, self.synthesis_config['letter_substitution_error_ratio'], len(oligo))
            oligo_with_letters_substitution = [''] * len(oligo)
            for letter_idx, letter in enumerate(oligo):
                if substitution_idx[letter_idx] == 1:
                    diff = {'A', 'C', 'G', 'T'} - set(letter)
                    oligo_with_letters_substitution[letter_idx] = random.choice(''.join(diff))
                else:
                    oligo_with_letters_substitution[letter_idx] = letter
            dna_list[row_idx] = ''.join(oligo_with_letters_substitution)
        return dna_list

    def get_x_list(self, payload: List[str]):
        x_list = []
        for z in payload:
            for key, val, in self.k_mer_representative_to_z.items():
                if val == z:
                    x_list.append(key)
        return x_list

    def constrained_sum_sample_pos(self, n, total):
        """Return a randomly chosen list of n positive integers summing to total.
        Each such list is equally likely to occur."""

        dividers = sorted(random.sample(range(1, total), n - 1))
        return [a - b for a, b in zip(dividers + [total], [0] + dividers)]
