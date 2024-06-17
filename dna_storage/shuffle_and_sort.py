import os
from pathlib import Path
import sqlite3

from unireedsolomon.unireedsolomon import RSCodecError

from dna_storage.config import PathLike
from dna_storage.rs_adapter import RSBarcodeAdapter


def shuffle(shuffle_db_file: PathLike, input_file: PathLike, output_file: PathLike):

    try:
        os.remove(shuffle_db_file)
    except OSError:
        pass

    if isinstance(shuffle_db_file, Path):
        shuffle_db_file = str(shuffle_db_file)
    conn = sqlite3.connect(shuffle_db_file)
    c = conn.cursor()

    c.execute('''CREATE TABLE shuffle_table
                 (idx real, data text)''')

    with open(input_file, 'r') as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            c.execute("INSERT INTO shuffle_table VALUES (" + str(idx) + ", '" + line + "')")

    conn.commit()

    c.execute("SELECT * FROM shuffle_table ORDER BY RANDOM()")
    with open(output_file, 'w+') as f:
        for idx, line in c:
            f.write(line + '\n')


def sample_oligos_from_file(input_file: PathLike, output_file: PathLike, number_of_oligos: int, number_of_blocks: int = 1):
    number_of_sampled_oligo = number_of_oligos * number_of_blocks
    with open(input_file, 'r') as input_file, open(output_file, 'w+') as output_file:
        for idx, line in enumerate(input_file):
            if idx < number_of_sampled_oligo:
                output_file.write(line)
            else:
                return


def sort_oligo_file(barcode_len: int, barcode_rs_len: int,
                    sort_db_file: PathLike, input_file: PathLike, output_file: PathLike,
                    barcode_coder: RSBarcodeAdapter):
    try:
        os.remove(sort_db_file)
    except OSError:
        pass

    if isinstance(sort_db_file, Path):
        sort_db_file = str(sort_db_file)
    conn = sqlite3.connect(sort_db_file)
    c = conn.cursor()

    c.execute('''CREATE TABLE sort_table
                     (barcode text, data text)''')

    with open(input_file, 'r') as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            barcode = line[:barcode_len+barcode_rs_len]
            payload = line[barcode_len+barcode_rs_len:]

            barcode_list = [c for c in barcode]
            try:
                barcode_decoded = barcode_coder.decode(barcode_encoded=barcode_list)
            except RSCodecError:
                continue

            barcode_decoded = ''.join(barcode_decoded)
            c.execute("INSERT INTO sort_table VALUES ('" + barcode_decoded + "', '" + payload + "')")

    c.execute("SELECT * FROM sort_table ORDER BY barcode")
    with open(output_file, 'w+') as f:
        for barcode, payload in c:
            f.write(barcode + payload + '\n')
