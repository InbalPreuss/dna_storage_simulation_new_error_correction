# Use sagemath and sage to use the bch code to encode and decode the data.
import math

import bchlib as bch

import itertools

import numpy as np
from unireedsolomon.unireedsolomon import rs, RSCodecError
from unireedsolomon.unireedsolomon import ff


class BCHPayloadAdapter:
    def __init__(self, bits_per_z, payload_len, payload_rs_len):
        self.bits_per_z = bits_per_z
        self.payload_len = payload_len
        # self.n = payload_len + payload_rs_len
        self.n = (payload_len + payload_rs_len)
        # self.n = (payload_len + payload_rs_len) * 3
        # self.n = 15
        # self.k = payload_len
        # self.k = payload_len * 3
        self.k = payload_len
        self.alphabet = list(range((payload_len + payload_rs_len)))  # Alphabet is {0,1,2,3,4,5,6,7}
        # self.c_exp = bits_per_z
        # self.c_exp = math.ceil(math.log2(self.n + 1))
        self.c_exp = math.ceil(math.log2(self.n))
        self.generator = 3
        self.t = (self.n - self.k) // 2

        # Initialize the BCH coder
        self.payload_coder = bch.BCH(self.t, m=self.c_exp)



        # self.bits_per_z = bits_per_z
        # self.payload_len = payload_len
        # self.n = payload_len + payload_rs_len
        # self.k = payload_len
        # self.alphabet = list(range(self.n))  # Alphabet is {0,1,2,3,4,5,6,7}
        # self.c_exp = bits_per_z
        # self.c_exp = 5
        # self.generator = 3
        # self.t = (self.n - self.k) // 2
        #
        # # Initialize the BCH coder
        # self.payload_coder = bch.BCH(self.t, m=self.c_exp)

        x=4



        # self.bits_per_z = bits_per_z
        # self.payload_len = payload_len
        # n = payload_len + payload_rs_len
        # k = payload_len
        # alphabet = list(range(n))
        # c_exp = bits_per_z
        # generator = 3
        # # Calculated error-correcting capability
        # t = (n - k) // 2
        # # t = (n - k)
        # # prim = ff.find_prime_polynomials(generator=generator, c_exp=c_exp, fast_primes=False, single=True)
        #
        # self.payload_coder = bch.BCH(t, m=c_exp)
    def encode(self, payload):
        # ff.set_globals(*self.ff_globals)
        # payload_as_int = [self._payload_to_int[z] for z in payload]
        # ff.set_globals(*self.ff_globals)
        # # Convert payload to bytes if necessary
        # if isinstance(payload, str):
        #     payload = payload.encode('utf-8')

        # Convert each digit to its 3-bit binary representation
        binary_representation = ''.join(format(int(char), '03b') for char in payload)

        # Function to convert binary string to bytearray
        def binary_string_to_bytearray(binary_string):
            # Split the binary string into chunks of 8 bits
            byte_chunks = [binary_string[i:i + 8] for i in range(0, len(binary_string), 8)]
            # Convert each chunk to a byte and create a bytearray
            return bytearray(int(chunk, 2) for chunk in byte_chunks)

        # Convert the binary representation to a bytearray
        ecc_byte_array = binary_string_to_bytearray(binary_representation)

        binary_string = ''.join(format(byte, '08b') for byte in ecc_byte_array)

        # Split the binary string into 3-bit groups
        three_bit_groups = [binary_string[i:i + 3] for i in range(0, len(binary_string), 3)]

        # Ensure the ECC is within the alphabet
        # payload_ecc_encoded = self._ensure_alphabet(payload_ecc_encoded)

        # Combine the payload and ECC
        # payload_encoded = payload + payload_ecc_encoded
        payload_encoded = payload + ecc_byte_array
        return payload_encoded



        # payload_ecc_encoded = self.payload_coder.encode(payload)
        # payload_encoded = payload + payload_ecc_encoded
        # return payload_encoded

    def decode(self, payload_encoded):
        ff.set_globals(*self.ff_globals)
        payload_as_int = [self._payload_to_int[z] for z in payload_encoded]
        if self._payload_coder.check_fast(payload_as_int):
            return payload_encoded[0:self.payload_len]
        else:
            try:
                payload_as_gf, rs_as_gf = self._payload_coder.decode(payload_as_int, nostrip=True, return_string=False)
            except RSCodecError:
                return payload_encoded[0:self.payload_len]
            payload = [self._int_to_payload[i] for i in payload_as_gf]
            return payload

#
# class RSPayloadAdapter:
#     def __init__(self, bits_per_z, payload_len, payload_rs_len):
#         ###########################
#         bits_per_z = 3
#         payload_len = 4
#         payload_rs_len = 4
#         ###########################
#         self.bits_per_z = bits_per_z
#         self.payload_len = payload_len
#         alphabet = ['Z{}'.format(i) for i in range(1, 2**bits_per_z+1)]
#         n = payload_len + payload_rs_len
#         k = payload_len
#         c_exp = bits_per_z
#         generator = 2
#         prim = ff.find_prime_polynomials(generator=generator, c_exp=c_exp, fast_primes=False, single=True)
#         # prim = ff.find_prime_polynomials(generator=generator, c_exp=c_exp, fast_primes=False, single=False)
#
#         self._payload_coder = rs.RSCoder(n=n, k=k, generator=generator, prim=prim, c_exp=c_exp)
#         self.ff_globals = ff.get_globals()
#         self._payload_to_int = {''.join(vv): i for i, vv in enumerate(alphabet)}
#         self._int_to_payload = {i: vv for vv, i in self._payload_to_int.items()}
#
#     def encode(self, payload):
#         ff.set_globals(*self.ff_globals)
#         payload_as_int = [self._payload_to_int[z] for z in payload]
#         ff.set_globals(*self.ff_globals)
#         payload_encoded_as_polynomial = self._payload_coder.encode_fast(payload_as_int, return_string=False)
#         payload_encoded = [self._int_to_payload[z] for z in payload_encoded_as_polynomial]
#         return payload_encoded
#
#     def decode(self, payload_encoded):
#         ff.set_globals(*self.ff_globals)
#         payload_as_int = [self._payload_to_int[z] for z in payload_encoded]
#         if self._payload_coder.check_fast(payload_as_int):
#             return payload_encoded[0:self.payload_len]
#         else:
#             try:
#                 payload_as_gf, rs_as_gf = self._payload_coder.decode(payload_as_int, nostrip=True, return_string=False)
#             except RSCodecError:
#                 return payload_encoded[0:self.payload_len]
#             payload = [self._int_to_payload[i] for i in payload_as_gf]
#             return payload


BCHWideAdapter = BCHPayloadAdapter


if __name__ == '__main__':
    # data = '1234'
    #
    # bits_per_z = 6
    # payload_len = 4
    # payload_rs_len = 4
    # bch_adapter = BCHWideAdapter(bits_per_z, payload_len, payload_rs_len)
    # payload_encoded = bch_adapter.encode(payload=data)
    # data_ecc = bytearray(payload_encoded)
    # # insert error to data_ecc
    # data_ecc[0] ^= 0x01
    # data = data_ecc[:-bch_adapter.payload_coder.ecc_bytes]
    # ecc_dec = data_ecc[-bch_adapter.payload_coder.ecc_bytes:]
    # nerr = bch_adapter.payload_coder.decode(data, ecc_dec)


    n = 6
    k = 2
    c_exp = 6 # bits_per_z
    generator = 2
    prim = 11

    # Calculated error-correcting capability
    t = (n - k) // 2

    # Initialize the BCH encoder/decoder with the given parameters
    bch = bch.BCH(t, m=c_exp)

    # bch = bchlib.BCH(16, m=13)
    data = b'1234567'
    ecc = bch.encode_message(data)
    data_ecc = data + ecc
    data_ecc = bytearray(data_ecc)
    #insert error to data_ecc
    data_ecc[0] ^= 0x01
    data = data_ecc[:-bch.ecc_bytes]
    ecc_dec = data_ecc[-bch.ecc_bytes:]
    nerr = bch.decode(data, ecc_dec)
    print(f'data={data}')
    print('nerr:', nerr)
    print('syn:', bch.syn)
    print('errloc:', bch.errloc)

    d = bch.correct(data, ecc_dec)
    print(f'data={data}')
