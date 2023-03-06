import os
from random import choice
from string import ascii_letters

import numpy as np

from dna_storage.config import PathLike


class TextFileToBinaryFile:
    def __init__(self, input_file: str,
                 output_file: str,
                 payload_len: int,
                 bits_per_z: int,
                 oligos_per_block_len: int,
                 k_mer: int):
        self.input_file = input_file
        self.output_file = output_file
        self.payload_len = payload_len
        self.bits_per_z = bits_per_z
        self.oligos_per_block_len = oligos_per_block_len
        self.k_mer = k_mer

    def run(self):
        with open(self.input_file, 'r', encoding='utf-8') as input_file, open(self.output_file, 'w', encoding='utf-8') as output_file:
            oligo_len_binary = int(self.payload_len * self.bits_per_z)
            accumulation = ''
            number_of_binary_oligos_written = 0
            for line in input_file:
                text_data = line
                binary_data = text_to_bits(text_data)
                accumulation += binary_data
                while len(accumulation) >= oligo_len_binary:
                    to_write = accumulation[:oligo_len_binary]
                    accumulation = accumulation[oligo_len_binary:]
                    output_file.write(to_write + '\n')
                    number_of_binary_oligos_written += 1
            z_fill = 0

            # pad the last oligo to have length "oligo_len_binary"
            if len(accumulation) > 0:
                binary_data_padded, z_fill = self.transform_text_to_binary_string(binary_data=accumulation)
                output_file.write(binary_data_padded + '\n')
                number_of_binary_oligos_written += 1

            # pad to a multiplication of "oligos_per_block_for_rs"
            zeros_block, number_of_missing_rows_to_block = self.zero_pad_to_blocks_of_size(number_of_binary_oligos_written=number_of_binary_oligos_written, oligo_len_binary=oligo_len_binary)
            if zeros_block != '':
                output_file.write(zeros_block + '\n')

            n_zeros = (number_of_missing_rows_to_block * oligo_len_binary) + z_fill
            z_fill_text = "{0:b}".format(n_zeros).rjust(oligo_len_binary, '0')
            output_file.write(z_fill_text + '\n')

    def transform_text_to_binary_string(self, binary_data: str):
        oligo_len_binary = int(self.payload_len * self.bits_per_z)
        binary_data_len = len(binary_data)

        number_of_binary_oligos = np.ceil(binary_data_len / oligo_len_binary)
        total_binary_len = int(number_of_binary_oligos * oligo_len_binary)

        binary_data_padded = binary_data.ljust(total_binary_len, '0')
        z_fill = total_binary_len - binary_data_len

        return binary_data_padded, z_fill

    def zero_pad_to_blocks_of_size(self, number_of_binary_oligos_written: int, oligo_len_binary: int):
        excess_lines = number_of_binary_oligos_written % self.oligos_per_block_len
        number_of_missing_rows_to_block = self.oligos_per_block_len - excess_lines - 1
        # -1 because we write an extra lines. the number of zeros we appended to the last line of real data
        zeros_array = ['0'*oligo_len_binary for _ in range(number_of_missing_rows_to_block)]
        zeros_block = '\n'.join(zeros_array)
        return zeros_block, number_of_missing_rows_to_block


class DecoderResultToBinary:
    def __init__(self, input_file: PathLike,
                 output_file: PathLike,
                 barcode_len: int) -> None:

        self.input_file = input_file
        self.output_file = output_file
        self.barcode_len = barcode_len

    def run(self) -> None:
        with open(self.input_file, 'r', encoding='utf-8') as input_file, open(self.output_file, 'w', encoding='utf-8') as output_file:
            for idx, line in enumerate(input_file):
                barcode_and_payload = line.strip()
                barcode, payload = barcode_and_payload[:self.barcode_len], barcode_and_payload[
                                                        self.barcode_len:]
                output_file.write(payload + '\n')


class BinaryResultToText:
    def __init__(self, input_file: PathLike,
                 output_file: PathLike,
                 barcode_len: int,
                 payload_len: int,
                 bits_per_z: int) -> None:

        self.input_file = input_file
        self.output_file = output_file
        self.barcode_len = barcode_len
        self.payload_len = payload_len
        self.bits_per_z = bits_per_z
        open(self.output_file, 'w').close()

    def run(self) -> None:
        oligo_len_binary = int(self.payload_len * self.bits_per_z)
        with open(self.input_file, 'r+', encoding='utf-8') as input_file:
            if os.name == 'nt':
                newline_size = 2
            else:
                newline_size = 1
            failed_on_seek = False
            try:
                input_file.seek(0, os.SEEK_END)
                input_file.seek(input_file.tell() - oligo_len_binary - 1*newline_size, os.SEEK_SET)
            except ValueError:
                input_file.seek(0)
                failed_on_seek = True

            if not failed_on_seek:
                for idx, line in enumerate(input_file):
                    if idx == 0:
                        payload = line.strip()
                        try:
                            z_fill = int(payload, 2)
                            padded_rows_with_zeros, padded_zeros_in_first_line = divmod(z_fill, oligo_len_binary)
                            last_line_backslash_n = 1
                            truncate_length = (padded_zeros_in_first_line + (last_line_backslash_n * newline_size) +
                                               (padded_rows_with_zeros + 1) * (oligo_len_binary + newline_size))
                            input_file.seek(0, os.SEEK_END)
                            input_file.seek(input_file.tell() - truncate_length, os.SEEK_SET)
                            input_file.truncate()
                        except ValueError:
                            pass
                        input_file.seek(0)
                        break

        with open(self.input_file, 'r+', encoding='utf-8') as input_file, open(self.output_file, 'w', encoding='utf-8') as output_file:
            accumulation = ''
            utf_chars_sizes = [32, 24, 16, 8]
            for idx, line in enumerate(input_file):
                payload = line.strip()

                accumulation += payload
                stop = False
                while len(accumulation) >= utf_chars_sizes[0] and stop is False:
                    for size in utf_chars_sizes:
                        try:
                            bits = accumulation[:size]
                            text = text_from_bits(bits)
                            accumulation = accumulation[size:]
                            if text != '\x00':
                                output_file.write(text)
                            break
                        except UnicodeDecodeError:
                            if size == utf_chars_sizes[-1]:
                                stop = True
                                break

            if len(accumulation) > 0:
                try:
                    text = text_from_bits(accumulation)
                    output_file.write(text)
                except UnicodeDecodeError:
                    pass


def generate_random_text_file(size_kb: int, file: PathLike) -> None:
    text = ''.join(choice(ascii_letters) for i in range(1024*size_kb))
    with open(file, 'w') as f:
        f.write(text)


def text_to_bits(text: str, encoding: str = 'utf-8', errors: str = 'surrogatepass') -> str:
    bits = bin(int.from_bytes(text.encode(encoding, errors), 'big'))[2:]
    return bits.zfill(8 * ((len(bits) + 7) // 8))


def text_from_bits(bits: str, encoding: str = 'utf-8', errors: str = 'surrogatepass') -> str:
    n = int(bits, 2)
    return n.to_bytes((n.bit_length() + 7) // 8, 'big').decode(encoding, errors) or '\0'
