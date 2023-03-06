
from dna_storage.reedsolomon import barcode_rs_encode, barcode_rs_decode
from dna_storage.reedsolomon import rs4096_encode, rs4096_decode


def test_reed_solomon_z_encode_decode():
    z_list = ['Z1' for i in range(120)]
    z_encoded = rs4096_encode(z_list)

    n_err = 7
    z_encoded_error = z_encoded[:120-n_err]
    z_encoded_error = ['Z2' for i in range(n_err)] + z_encoded_error
    z_encoded_error.extend(z_encoded[120:])
    z_list_decoded = rs4096_decode(z_encoded_error, verify_only=False)

    assert ''.join(z_list_decoded) == ''.join(z_list)

    n_err = 8
    z_encoded_error = z_encoded[:120 - n_err]
    z_encoded_error = ['Z2' for i in range(n_err)] + z_encoded_error
    z_encoded_error.extend(z_encoded[120:])
    z_list_decoded = rs4096_decode(z_encoded_error, verify_only=False)

    assert ''.join(z_list_decoded) != ''.join(z_list)


def test_reed_solomon_barcode_encode_decode():
    barcode = 'AAAAAAAAAAAA'
    barcode_list = [i for i in barcode]
    barcode_encoded = barcode_rs_encode(barcode)

    barcode_encoded_error = [i for i in 'CCAAAAAAAAAA' + ''.join(barcode_encoded[12:])]
    barcode_decoded = barcode_rs_decode(barcode_encoded_error, verify_only=False)
    assert ''.join(barcode_list) == ''.join(barcode_decoded)

    barcode_encoded_error_wrong = [i for i in 'CACAAAAAAAAA' + ''.join(barcode_encoded[12:])]
    barcode_decoded_wrong = barcode_rs_decode(barcode_encoded_error_wrong, verify_only=False)
    assert ''.join(barcode_list) != ''.join(barcode_decoded_wrong)


if __name__ == '__main__':
    # test_reed_solomon_z_encode_decode()
    test_reed_solomon_barcode_encode_decode()
