import itertools
from collections import Counter
import re
from typing import Union, Dict, List
from pathlib import Path

from unireedsolomon import RSCodecError

from dna_storage.rs_adapter import RSBarcodeAdapter, RSPayloadAdapter, RSWideAdapter
from dna_storage import utils
from dna_storage.utils import chunker


#################################################################
# @ Class: Decoder
# @ Description: Retrieve the oligo to the oligo that was written
#                in originally
#################################################################
class Decoder:
    def __init__(self, barcode_len: int,
                 barcode_total_len: int,
                 payload_len: int,
                 payload_total_len: int,
                 input_file: str,
                 shrink_dict: Dict,
                 min_number_of_oligos_per_barcode: int,
                 k_mer: int,
                 k_mer_representative_to_z: Dict,
                 z_to_binary: Dict,
                 subset_size: int,
                 oligos_per_block_len: int,
                 oligos_per_block_rs_len: int,
                 drop_if_not_exact_number_of_chunks: bool,
                 barcode_coder: RSBarcodeAdapter,
                 payload_coder: RSPayloadAdapter,
                 wide_coder: RSWideAdapter,
                 results_file: Union[Path, str],
                 results_file_z_after_rs_wide: Union[Path, str],
                 results_file_z_before_rs_payload: Union[Path, str],
                 results_file_z_after_rs_payload: Union[Path, str],
                 ):
        self.input_file = input_file
        self.barcode_len = barcode_len
        self.barcode_total_len = barcode_total_len
        self.payload_len = payload_len
        self.payload_total_len = payload_total_len
        self.shrink_dict = shrink_dict
        self.min_number_of_oligos_per_barcode = min_number_of_oligos_per_barcode
        self.k_mer = k_mer
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.z_to_binary = z_to_binary
        self.subset_size = subset_size
        self.oligos_per_block_len = oligos_per_block_len
        self.oligos_per_block_rs_len = oligos_per_block_rs_len
        self.drop_if_not_exact_number_of_chunks = drop_if_not_exact_number_of_chunks
        self.results_file = results_file
        self.results_file_z_before_rs_payload = results_file_z_before_rs_payload
        self.results_file_z_after_rs_payload = results_file_z_after_rs_payload
        self.results_file_z_after_rs_wide = results_file_z_after_rs_wide
        open(self.results_file, 'w').close()
        open(self.results_file_z_before_rs_payload, 'w').close()
        open(self.results_file_z_after_rs_payload, 'w').close()
        open(self.results_file_z_after_rs_wide, 'w').close()
        self.barcode_generator = utils.dna_sequence_generator(sequence_len=self.barcode_len)
        self.barcode_coder = barcode_coder
        self.payload_coder = payload_coder
        self.wide_coder = wide_coder
        self.payload_total_len_nuc = (payload_total_len * k_mer)

    def run(self):
        barcode_prev = ''
        payload_accumulation = []
        dummy_payload = ['Z1' for _ in range(self.payload_len)]
        total_oligos_per_block_with_rs_oligos = self.oligos_per_block_len + self.oligos_per_block_rs_len
        with open(self.input_file, 'r', encoding='utf-8') as file:
            unique_payload_block_with_rs = []
            unique_barcode_block_with_rs = []
            for idx, line in enumerate(file):
                barcode_and_payload = line.split(sep=' ')[0].rstrip()
                barcode = barcode_and_payload[:self.barcode_len]
                payload = barcode_and_payload[self.barcode_len:]

                if barcode != barcode_prev:
                    next_barcode_should_be = "".join(next(self.barcode_generator))
                    if next_barcode_should_be != barcode:
                        unique_payload_block_with_rs.append(dummy_payload)
                        unique_barcode_block_with_rs.append(next_barcode_should_be)
                    if len(payload_accumulation) > self.min_number_of_oligos_per_barcode:
                        unique_payload = self.dna_to_unique_payload(payload_accumulation=payload_accumulation)
                        self.save_z_before_rs(barcode=barcode_prev, payload=unique_payload)
                        unique_payload_corrected = self.error_correction_payload(payload=unique_payload)
                        self.save_z_after_rs(barcode=barcode_prev, payload=unique_payload_corrected)
                        if len(unique_payload_corrected) > 0:
                            unique_payload_block_with_rs.append(unique_payload_corrected)
                        else:
                            unique_payload_block_with_rs.append(dummy_payload)

                        unique_barcode_block_with_rs.append(barcode_prev)
                        if len(unique_payload_block_with_rs) >= total_oligos_per_block_with_rs_oligos:
                            self.save_block_to_binary(
                                unique_barcode_block_with_rs[:total_oligos_per_block_with_rs_oligos],
                                unique_payload_block_with_rs[:total_oligos_per_block_with_rs_oligos]
                            )
                            unique_payload_block_with_rs = unique_payload_block_with_rs[total_oligos_per_block_with_rs_oligos:]
                            unique_barcode_block_with_rs = unique_barcode_block_with_rs[total_oligos_per_block_with_rs_oligos:]
                    payload_accumulation = [payload]
                    barcode_prev = barcode
                else:
                    payload_accumulation.append(payload)

            if len(payload_accumulation) > self.min_number_of_oligos_per_barcode:
                unique_payload = self.dna_to_unique_payload(payload_accumulation=payload_accumulation)
                self.save_z_before_rs(barcode=barcode_prev, payload=unique_payload)
                unique_payload_corrected = self.error_correction_payload(payload=unique_payload)
                self.save_z_after_rs(barcode=barcode_prev, payload=unique_payload_corrected)
                if len(unique_payload_corrected) > 0:
                    unique_payload_block_with_rs.append(unique_payload_corrected)
                else:
                    unique_payload_block_with_rs.append(dummy_payload)
                unique_barcode_block_with_rs.append(barcode_prev)
                while len(unique_payload_block_with_rs) < total_oligos_per_block_with_rs_oligos:
                    unique_payload_block_with_rs.append(dummy_payload)
                    next_barcode_should_be = "".join(next(self.barcode_generator))
                    unique_barcode_block_with_rs.append(next_barcode_should_be)
                self.save_block_to_binary(
                    unique_barcode_block_with_rs[:total_oligos_per_block_with_rs_oligos],
                    unique_payload_block_with_rs[:total_oligos_per_block_with_rs_oligos]
                )

    def dna_to_unique_payload(self, payload_accumulation: List[str]) -> List[str]:
        shrunk_payload = self.shrink_payload(payload_accumulation=payload_accumulation)
        shrunk_payload_histogram = self.payload_histogram(payload=shrunk_payload)
        unique_payload = self.payload_histogram_to_payload(payload_histogram=shrunk_payload_histogram)
        return unique_payload

    def save_block_to_binary(self, unique_barcode_block_with_rs: List[str],
                             unique_payload_block_with_rs: List[List[str]]) -> None:
        unique_payload_block = self.wide_rs(unique_payload_block_with_rs)
        for unique_barcode, unique_payload in zip(unique_barcode_block_with_rs, unique_payload_block):
            binary = self.unique_payload_to_binary(payload=unique_payload)
            self.save_z_after_rs_wide(barcode=unique_barcode, payload=unique_payload)
            self.save_binary(binary=binary, barcode_prev=unique_barcode)

    def wide_rs(self, unique_payload_block_with_rs):
        rs_removed = [[] for _ in range(int(self.oligos_per_block_len))]
        for col in range(len(unique_payload_block_with_rs[0])):
            payload = [elem[col] for elem in unique_payload_block_with_rs]
            col_without_rs = self.error_correction_payload(payload=payload, payload_or_wide='wide')
            if len(col_without_rs) > self.oligos_per_block_len:
                import logging
                logger = logging.getLogger()
                logger.error(f"DECODER ERROR, "
                             f"unique_payload_block_with_rs: {unique_payload_block_with_rs}, "
                             f"payload: {payload}, "
                             f"col_without_rs: {col_without_rs}")
            for idx, z in enumerate(col_without_rs):
                rs_removed[idx].append(z)
        return rs_removed

    def unique_payload_to_binary(self, payload: List[str]) -> str:
        binary = []
        for z in payload:
            try:
                binary.append(self.z_to_binary[z])
            except KeyError:
                return ''
        binary = ["".join(map(str, tup)) for tup in binary]
        return "".join(binary)

    def wrong_barcode_and_payload_len(self, barcode_and_payload: str) -> bool:
        return len(barcode_and_payload) != self.barcode_len + self.payload_total_len_nuc

    def shrink_payload(self, payload_accumulation: List[str]) -> List[List[str]]:
        if self.k_mer == 1:
            return [payload_accumulation]
        k_mer_accumulation = []
        # When we inserted errors per nuc (not per entire Z). we used those 3 algorithms to fix some of the errors.
        # for payload in payload_accumulation:
        #     k_mer_list = self.get_transformed_oligo_with_correct_len(payload)
        #     k_mer_accumulation.append(k_mer_list)
        # return k_mer_accumulation
        for payload in payload_accumulation:
            k_mer_list = []
            oligo_valid = True
            if self.drop_if_not_exact_number_of_chunks:
                if not self.payload_total_len_nuc - (self.k_mer - 1) <= len(payload) <= self.payload_total_len_nuc + (self.k_mer - 1):
                    continue
            else:
                payload = payload.ljust(self.payload_total_len_nuc, 'R')[:self.payload_total_len_nuc]
            for k_letters in chunker(payload, size=self.k_mer):
                k_mer_list.append(self.shrink_dict.get(k_letters, "Xdummy"))
            if oligo_valid:
                k_mer_accumulation.append(k_mer_list)
        return k_mer_accumulation

    def get_transformed_oligo_with_correct_len(self, payload: str) -> List[str]:
        k_mer_list = []
        payload_len = len(payload)
        delta = payload_len - self.payload_total_len_nuc
        i = 0
        while True:
            k_letters = payload[i:i + self.k_mer]
            try:
                k_mer = self.shrink_dict[k_letters]
                k_mer_list.append(k_mer)
            except KeyError:
                if delta < 0:
                    delta += 1
                    i -= 1
                    k_mer_list.append("Xdummy")
                elif delta > 0:
                    delta -= 1
                    i += 1
                    try:
                        next_letter = payload[i + self.k_mer]
                    except IndexError:
                        k_mer_list.append("Xdummy")
                    else:
                        success = False
                        for letters in itertools.combinations(k_letters, self.k_mer - 1):
                            k_letters_try = "".join(letters) + next_letter
                            try:
                                k_mer = self.shrink_dict[k_letters_try]
                                k_mer_list.append(k_mer)
                                success = True
                            except KeyError:
                                pass
                        if not success:
                            k_mer_list.append("Xdummy")
                else:
                    k_mer_list.append("Xdummy")
            i += 3
            if len(k_mer_list) >= self.payload_total_len:
                return k_mer_list

    def payload_histogram(self, payload: List[List[str]]) -> List[Counter]:
        hist = []
        for col_idx in range(self.payload_total_len):
            col = [letter[col_idx] for letter in payload]
            letter_counts = Counter(col)
            hist.append(letter_counts)
        return hist

    def error_correction_payload(self, payload: Union[str, List[str]], payload_or_wide: str = 'payload') -> List[str]:
        if isinstance(payload, str):
            payload = [c for c in payload]

        if payload_or_wide == 'payload':
            payload_decoded = self.payload_coder.decode(payload_encoded=payload)
        else:
            payload_decoded = self.wide_coder.decode(payload_encoded=payload)

        return payload_decoded

    def error_correction_barcode(self, barcode: Union[str, List[str]]) -> str:
        if isinstance(barcode, str):
            barcode = [c for c in barcode]

        try:
            barcode_decoded = self.barcode_coder.decode(barcode_encoded=barcode)
        except RSCodecError:
            return ''

        if isinstance(barcode_decoded, list):
            barcode_decoded = ''.join(barcode_decoded)
        return barcode_decoded

    def payload_histogram_to_payload(self, payload_histogram: List[Counter]) -> List[str]:
        result_payload = []
        for counter in payload_histogram:
            del counter['Xdummy']
            reps = counter.most_common(self.subset_size)
            if len(reps) != self.subset_size:
                # BAD
                # return [0]
                # OK
                # reps = [('X' + str(i), i) for i in range(1, self.subset_size + 1)]
                # BEST
                s = set(['X' + str(i) for i in range(1, self.subset_size + 1)])
                t = set([r[0] for r in reps])
                diff = sorted(list(s-t))
                for missing_rep_idx in range(self.subset_size - len(reps)):
                    reps.append((diff[missing_rep_idx], 1))
            k_mer_rep = tuple(self.sorted_human([rep[0] for rep in reps]))
            try:
                z = self.k_mer_representative_to_z[k_mer_rep]
            except KeyError:
                z = 'Z1'  # The X tuple is out of range
            result_payload.append(z)

        return result_payload

    def save_z_before_rs(self, payload: List[str], barcode: str) -> None:
        with open(self.results_file_z_before_rs_payload, 'a+', encoding='utf-8') as f:
            f.write(barcode + "," + ",".join(payload) + '\n')

    def save_z_after_rs(self, payload: List[str], barcode: str) -> None:
        with open(self.results_file_z_after_rs_payload, 'a+', encoding='utf-8') as f:
            f.write(barcode + "," + ",".join(payload) + '\n')

    def save_z_after_rs_wide(self, payload: List[str], barcode: str) -> None:
        with open(self.results_file_z_after_rs_wide, 'a+', encoding='utf-8') as f:
            f.write(barcode + "," + ",".join(payload) + '\n')

    def save_binary(self, binary: str, barcode_prev: str) -> None:
        with open(self.results_file, 'a+', encoding='utf-8') as f:
            f.write(barcode_prev + binary + '\n')

    @staticmethod
    def sorted_human(iterable: List[str]) -> List[str]:
        """ Sort the given iterable in the way that humans expect."""
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(iterable, key=alphanum_key)
