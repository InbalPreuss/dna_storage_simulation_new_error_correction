import json
from pathlib import Path

from dna_storage.config import build_config

from dna_storage.decoder import Decoder
from dna_storage.shuffle_and_sort import sort_oligo_file
from dna_storage.text_handling import BinaryResultToText, DecoderResultToBinary

# binary_results_to_text = BinaryResultToText(
#     input_file="/Users/omri/git/DnaStorage/data/testing goood/check_final_6/simulation_data.11.binary_results_file.dna",
#                                                     output_file="a",
#                                                     barcode_len=12,
#                                                     payload_len=120,
#                                                     bits_per_z=9)
# binary_results_to_text.run()
from tests.distributed import build_runs


def rebuild_config_from_outputs():
    output_dir = Path(
        "data/testing goood/check/[ subset size 5, bits per z 12 ][ number of oligos per barcode   1000 ][ number of oligos sampled after synthesis   1000 ][ errors, substitution 0     , deletion 0     , insertion 0.01   ] trial  9"
    )
    json_file = [file for file in output_dir.iterdir() if file.suffix == ".json"][0]
    input_file = [file for file in output_dir.iterdir() if file.suffix == ".txt"][0]
    with open(json_file, "r") as f:
        config_for_run = json.load(f)
    config = build_config(
        subset_size=config_for_run["size"],
        bits_per_z=config_for_run["bits_per_z"],
        letter_substitution_error_ratio=config_for_run["substitution_error"],
        letter_deletion_error_ratio=config_for_run["deletion_error"],
        letter_insertion_error_ratio=config_for_run["insertion_error"],
        number_of_oligos_per_barcode=config_for_run["number_of_oligos_per_barcode"],
        number_of_sampled_oligos_from_file=config_for_run["number_of_sampled_oligos_from_file"],
        output_dir=output_dir,
        input_text_file=input_file,
    )
    return config


config = rebuild_config_from_outputs()

shrink_dict = config['shrink_dict']
decoder = Decoder(barcode_len=config['barcode_len'],
                  barcode_total_len=config['barcode_total_len'],
                  payload_len=config['payload_len'],
                  payload_total_len=config['payload_total_len'],
                  input_file=config['sort_oligo_results_file'],
                  shrink_dict=shrink_dict,
                  min_number_of_oligos_per_barcode=config['min_number_of_oligos_per_barcode'],
                  k_mer=config['k_mer'],
                  k_mer_representative_to_z=config['algorithm_config']['k_mer_representative_to_z'],
                  z_to_binary=config['algorithm_config']['z_to_binary'],
                  subset_size=config['algorithm_config']['subset_size'],
                  oligos_per_block_len=config['oligos_per_block_len'],
                  oligos_per_block_rs_len=config['oligos_per_block_rs_len'],
                  barcode_coder=config['barcode_coder'],
                  payload_coder=config['payload_coder'],
                  wide_coder=config['wide_coder'],
                  results_file=config['decoder_results_file'],
                  results_file_z_before_rs_payload=config['decoder_results_file_z_before_rs_payload'],
                  results_file_z_after_rs_payload=config['decoder_results_file_z_after_rs_payload'],
                  results_file_z_after_rs_wide=config['decoder_results_file_z_after_rs_wide'],
                  )
decoder.run()

decoder_results_to_binary = DecoderResultToBinary(input_file=config['decoder_results_file'],
                                                  output_file=config['binary_results_file'],
                                                  barcode_len=config['barcode_len'])
decoder_results_to_binary.run()

binary_results_to_text = BinaryResultToText(input_file=config['binary_results_file'],
                                            output_file=config['text_results_file'],
                                            barcode_len=config['barcode_len'],
                                            payload_len=config['payload_len'],
                                            bits_per_z=config['algorithm_config']['bits_per_z'])
binary_results_to_text.run()
