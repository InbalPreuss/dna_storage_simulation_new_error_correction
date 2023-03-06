import glob
import os
import pathlib

from Bio.SeqIO.QualityIO import FastqGeneralIterator


#################################################################
# @ Function: get_seq_id_offset
# @ Input: seq_id
# @ Description: Calculate the sequence id offset
# @ Return: offset
#################################################################
def get_seq_id_offset(seq_id):
    total_seq_index = seq_id
    multiply_offset = 9
    number_of_id_digits = 1
    offset = 0

    while (total_seq_index - multiply_offset) > 0:
        offset += (multiply_offset * number_of_id_digits)
        total_seq_index -= multiply_offset
        multiply_offset *= 10
        number_of_id_digits += 1

    offset += ((total_seq_index - 1) * number_of_id_digits)
    return offset


#################################################################
# @ class: FastqHandling
# @ Description: FastqHandling handles the .fastq file and makes
#                a .txt file with the sorted Oligos by the barcode
#################################################################
class FastqHandling:

    #################################################################
    # @ Function: __init__
    # @ Input: barcode_len, payload_len, file_name
    # @ Description: init class FastqHandling
    #################################################################
    def __init__(self, barcode_len: int, payload_len: int, file_name: str):
        fastq_input_file_name = pathlib.Path(r'data/input/' + file_name + '.fastq')

        if fastq_input_file_name.exists() == 0:
            raise NameError('The file name ' + file_name + ' does not exist')

        self.file_name = file_name
        self.file_full_name = fastq_input_file_name
        self.file_extension = fastq_input_file_name.suffix
        self.file_full_name_set_ids_output = "data/fastq_output/" + file_name + "_set_ids.txt"
        self.file_full_name_sorted_output = "data/fastq_output/" + file_name + "_sorted.dna"
        self.barcode_len = barcode_len
        self.payload_len = payload_len

    #################################################################
    # @ Function: set_oligo_id
    # @ Description: Set new ids to the Oligos,
    #                from 1,2,3,4.... to the amount of oligos in file
    #                and put the ids + Oligo into a .txt file
    #################################################################
    def set_oligo_id(self):
        seq_id = 1
        try:
            os.makedirs(os.path.dirname(self.file_full_name_set_ids_output), exist_ok=True)
            file_name_set_ids = open(self.file_full_name_set_ids_output, "w")
            for (title, sequence, quality) in FastqGeneralIterator(self.file_full_name):
                seq_and_id = str(seq_id) + " " + str(sequence) + "\n"
                file_name_set_ids.write(seq_and_id)
                seq_id += 1
        except FileExistsError:
            print('There is a problem opening the file ' + file_name_set_ids)

        file_name_set_ids.close()

    #################################################################
    # @ Function: sort_oligo
    # @ Description: Sort the Oligos by barcode and insert
    #                all the Oligos into a .txt file
    #################################################################
    def sort_oligo(self):
        try:
            os.makedirs(os.path.dirname(self.file_full_name_sorted_output), exist_ok=True)
            file_name_set_ids = open(self.file_full_name_set_ids_output, "r")
            file_name_sorted = open(self.file_full_name_sorted_output, "w")
        except FileExistsError:
            print('There is a problem opening the file ' + file_name_set_ids + " or " + file_name_sorted)

        # Insert every Oligo sequence, their id and barcode to a list
        id_by_barcode_list = []
        for line in file_name_set_ids:
            # line = id, seq
            id_and_seq = line.split()
            id_by_barcode_list.append([id_and_seq[0], id_and_seq[1][:self.barcode_len]])

        # Sort ids by Oligo barcode
        sorted_oligo_by_barcode = sorted(id_by_barcode_list, key=lambda x: x[1])

        # Sort the oligo sequence by barcode and insert to .txt file
        # We do not want to read every time the entire lines to get to the specific line that we want,
        # therefore we calculate the location of the line and jump to that location
        for seq_id, barcode in sorted_oligo_by_barcode:
            # line_offset = prev_line_number * (payload_len + space + new_line) + seq_id_offset + prev_line_number
            line_offset = (int(seq_id) - 1) * (self.payload_len + 1 + 1) + get_seq_id_offset(int(seq_id)) + (
                    int(seq_id) - 1)

            # Jumps to the Oligo line and writes the line to the file
            file_name_set_ids.seek(line_offset)
            line = file_name_set_ids.readline()
            file_name_sorted.write(line)

        file_name_set_ids.close()
        file_name_sorted.close()

    def parse_fastq(self):
        self.set_oligo_id()
        self.sort_oligo()

        print("Finished parsing the Fastq file")
        return self.file_full_name_sorted_output
