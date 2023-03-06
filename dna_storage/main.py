from dna_storage.fastq_handling import FastqHandling
from dna_storage.text_handling import TextFileToBinaryFile, DecoderResultToBinary, BinaryResultToText
from dna_storage.decoder import Decoder
from dna_storage.encoder import Encoder
from dna_storage.mock_synthesizer import Synthesizer
from dna_storage.shuffle_and_sort import shuffle, sort_oligo_file, sample_oligos_from_file


def main(config):
    if config['write_text_to_binary']:
        print(f"1. write_text_to_binary")
        text_file_to_binary = TextFileToBinaryFile(input_file=config['input_text_file'],
                                                   output_file=config['binary_file_name'],
                                                   payload_len=config['payload_len'],
                                                   bits_per_z=config['algorithm_config']['bits_per_z'],
                                                   oligos_per_block_len=config['oligos_per_block_len'],
                                                   k_mer=config['k_mer'])
        text_file_to_binary.run()

    # Encode
    if config['do_encode']:
        print(f"2. encode")
        shrink_dict = config['shrink_dict']
        encoder = Encoder(barcode_len=config['barcode_len'],
                          barcode_rs_len=config['barcode_rs_len'],
                          payload_len=config['payload_len'],
                          payload_rs_len=config['payload_rs_len'],
                          binary_file_name=config['binary_file_name'],
                          shrink_dict=shrink_dict,
                          k_mer=config['k_mer'],
                          k_mer_representative_to_z=config['algorithm_config']['k_mer_representative_to_z'],
                          binary_to_z=config['algorithm_config']['binary_to_z'],
                          subset_size=config['algorithm_config']['subset_size'],
                          oligos_per_block_len=config['oligos_per_block_len'],
                          oligos_per_block_rs_len=config['oligos_per_block_rs_len'],
                          bits_per_z=config['algorithm_config']['bits_per_z'],
                          barcode_coder=config['barcode_coder'],
                          payload_coder=config['payload_coder'],
                          wide_coder=config['wide_coder'],
                          results_file=config['encoder_results_file'],
                          results_file_without_rs_wide=config['encoder_results_file_without_rs_wide'])
        number_of_blocks = encoder.run()

    # Synthesize
    if config['do_synthesize']:
        print(f"3. synthesize")
        synthesizer = Synthesizer(input_file=config['encoder_results_file'],
                                  results_file=config['synthesis_results_file'],
                                  synthesis_config=config['synthesis'],
                                  barcode_total_len=config['barcode_total_len'],
                                  subset_size=config['algorithm_config']['subset_size'],
                                  k_mer_representative_to_z=config['algorithm_config']['k_mer_representative_to_z'],
                                  k_mer_to_dna=config['algorithm_config']['k_mer_to_dna'],
                                  k_mer=config['k_mer'],
                                  mode=config['mode'])
        synthesizer.synthesize()

    # Shuffling the sorted synthesis results
    if config['do_shuffle']:
        print(f"4. shuffle")
        shuffle(shuffle_db_file=config['shuffle_db_file'],
                input_file=config['synthesis_results_file'],
                output_file=config['shuffle_results_file'])

    # Sample from the shuffled synthesis results
    if config['do_sample_oligos_from_file']:
        print(f"5. sample oligos from file")
        sample_oligos_from_file(input_file=config['shuffle_results_file'],
                                output_file=config['sample_oligos_results_file'],
                                number_of_oligos=config['number_of_sampled_oligos_from_file'],
                                number_of_blocks=number_of_blocks)

    # Sorting the shuffled synthesis results
    if config['do_sort_oligo_file']:
        print(f"6. sort oligo file")
        sort_oligo_file(barcode_len=config['barcode_len'],
                        barcode_rs_len=config['barcode_rs_len'],
                        sort_db_file=config['sort_oligo_db_file'],
                        input_file=config['sample_oligos_results_file'],
                        output_file=config['sort_oligo_results_file'],
                        barcode_coder=config['barcode_coder'])

    # Parsing Fastq data
    if config['do_fastq_handling']:
        print(f"7. fastq handling")
        file_name_sorted = FastqHandling(barcode_len=config['barcode_len'],
                                         payload_len=config['payload_len'],
                                         file_name=config['fastq_file_name']).parse_fastq()
    # Decode
    if config['do_decode']:
        print(f"8. decode")
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
                          drop_if_not_exact_number_of_chunks=config['drop_if_not_exact_number_of_chunks'],
                          barcode_coder=config['barcode_coder'],
                          payload_coder=config['payload_coder'],
                          wide_coder=config['wide_coder'],
                          results_file=config['decoder_results_file'],
                          results_file_z_before_rs_payload=config['decoder_results_file_z_before_rs_payload'],
                          results_file_z_after_rs_payload=config['decoder_results_file_z_after_rs_payload'],
                          results_file_z_after_rs_wide=config['decoder_results_file_z_after_rs_wide'],
                          )
        decoder.run()

    if config['decoder_results_to_binary']:
        print(f"9. results to binary")
        decoder_results_to_binary = DecoderResultToBinary(input_file=config['decoder_results_file'],
                                                          output_file=config['binary_results_file'],
                                                          barcode_len=config['barcode_len'])
        decoder_results_to_binary.run()

    if config['binary_results_to_text']:
        print(f"10. binary results to text")
        binary_results_to_text = BinaryResultToText(input_file=config['binary_results_file'],
                                                    output_file=config['text_results_file'],
                                                    barcode_len=config['barcode_len'],
                                                    payload_len=config['payload_len'],
                                                    bits_per_z=config['algorithm_config']['bits_per_z'])
        binary_results_to_text.run()


if __name__ == "__main__":
    from dna_storage.config import build_config
    config = build_config()
    main(config=config)
