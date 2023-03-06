from textwrap import wrap
from typing import Union, Dict, List, Tuple
from pathlib import Path

from dna_storage.rs_adapter import RSBarcodeAdapter, RSPayloadAdapter, RSWideAdapter
from dna_storage import utils


#################################################################
# @ Class: Encoder
# @ Description: Retrieve the oligo to the oligo that was written
#                in originally
#################################################################

class Encoder:
    def __init__(self, barcode_len: int,
                 barcode_rs_len: int,
                 payload_len: int,
                 payload_rs_len: int,
                 binary_file_name: str,
                 shrink_dict: Dict,
                 k_mer: int,
                 k_mer_representative_to_z: Dict,
                 binary_to_z: Dict,
                 subset_size: int,
                 bits_per_z: int,
                 oligos_per_block_len: int,
                 oligos_per_block_rs_len: int,
                 barcode_coder: RSBarcodeAdapter,
                 payload_coder: RSPayloadAdapter,
                 wide_coder: RSWideAdapter,
                 results_file: Union[Path, str],
                 results_file_without_rs_wide: Union[Path, str],
                 ):
        self.file_name = binary_file_name
        self.barcode_len = barcode_len
        self.barcode_rs_len = barcode_rs_len
        self.payload_len = payload_len
        self.payload_rs_len = payload_rs_len
        self.shrink_dict = shrink_dict
        self.k_mer = k_mer
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.binary_to_z_dict = binary_to_z
        self.subset_size = subset_size
        self.bits_per_z = bits_per_z
        self.oligos_per_block_len = oligos_per_block_len
        self.oligos_per_block_rs_len = oligos_per_block_rs_len
        self.results_file = results_file
        self.results_file_without_rs_wide = results_file_without_rs_wide
        open(self.results_file, 'w').close()
        self.barcode_generator = utils.dna_sequence_generator(sequence_len=self.barcode_len)
        self.barcode_coder = barcode_coder
        self.payload_coder = payload_coder
        self.wide_coder = wide_coder

    def run(self):
        number_of_blocks = 0
        with open(self.file_name, 'r', encoding='utf-8') as file:
            z_list_accumulation_per_block = []
            for line in file:
                line = line.strip('\n')
                z_list = []
                for binary_to_transform in wrap(line, self.bits_per_z):
                    z = self.binary_to_z(binary=binary_to_transform)
                    z_list.append(z)
                z_list_accumulation_per_block.append(z_list)
                if len(z_list_accumulation_per_block) == self.oligos_per_block_len:
                    number_of_blocks += 1
                    z_list_accumulation_with_rs = self.wide_block_rs(z_list_accumulation_per_block)
                    amount_oligos_per_block_len_to_write = self.oligos_per_block_len
                    for z_list in z_list_accumulation_with_rs:
                        oligo = self.z_to_oligo(z_list)
                        self.save_oligo(results_file=self.results_file, oligo=oligo)
                        if amount_oligos_per_block_len_to_write > 0:
                            self.save_oligo(results_file=self.results_file_without_rs_wide, oligo=oligo)
                            amount_oligos_per_block_len_to_write = amount_oligos_per_block_len_to_write - 1
                    z_list_accumulation_per_block = []
        return number_of_blocks

    def binary_to_z(self, binary: str) -> str:
        binary_tuple = tuple([int(b) for b in binary])
        return self.binary_to_z_dict[binary_tuple]

    def z_to_oligo(self, z_list: List[str]) -> str:
        oligo = self.add_payload_rs_symbols_for_error_correction(payload=z_list)
        barcode = next(self.barcode_generator)
        barcode = self.add_barcode_rs_symbols_for_error_correction(barcode=barcode)
        barcode = "".join(barcode)
        oligo.insert(0, barcode)
        return ",".join(oligo)

    def wide_block_rs(self, z_list_accumulation_per_block: List[List[str]]) -> List[List[str]]:
        rs_append = [[] for _ in range(int(self.oligos_per_block_len + self.oligos_per_block_rs_len))]
        for col in range(len(z_list_accumulation_per_block[0])):
            z_list = [elem[col] for elem in z_list_accumulation_per_block]
            col_with_rs = self.add_payload_rs_symbols_for_error_correction(payload=z_list, payload_or_wide='wide')
            for idx, z in enumerate(col_with_rs):
                rs_append[idx].append(z)
        return rs_append

    def add_payload_rs_symbols_for_error_correction(self, payload: Union[str, List[str]],
                                                    payload_or_wide: str = 'payload') -> List[str]:
        if isinstance(payload, str):
            payload = [c for c in payload]
        if payload_or_wide == 'payload':
            payload_encoded = self.payload_coder.encode(payload)
        else:
            payload_encoded = self.wide_coder.encode(payload)

        return payload_encoded

    def add_barcode_rs_symbols_for_error_correction(self, barcode: Tuple[str]) -> List[str]:
        barcode = list(barcode)
        barcode_encoded = self.barcode_coder.encode(barcode=barcode)
        return barcode_encoded

    def save_oligo(self, results_file: Union[Path, str], oligo: str) -> None:
        with open(results_file, 'a+', encoding='utf-8') as f:
            f.write(oligo + '\n')
