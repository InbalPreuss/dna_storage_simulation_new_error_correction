import itertools
from typing import Dict

import numpy as np
from unireedsolomon.unireedsolomon import rs, RSCodecError
from unireedsolomon.unireedsolomon import ff
from dna_storage.vt_syndrome import VTSyndrome
from dna_storage import utils as uts

class RSBarcodeAdapter:
    def __init__(self, bits_per_z, barcode_len, barcode_rs_len):
        self.bits_per_z = bits_per_z
        self._barcode_len = barcode_len
        alphabet = list(itertools.product('ACGT', 'ACGT'))
        n = int((barcode_len + barcode_rs_len) / 2)
        k = int(barcode_len / 2)
        c_exp = int(np.log2(len(alphabet)))
        generator = 3
        prim = ff.find_prime_polynomials(generator=generator, c_exp=c_exp, fast_primes=False, single=True)

        self._barcode_coder = rs.RSCoder(n=n, k=k, generator=generator, prim=prim, c_exp=c_exp)
        self.ff_globals = ff.get_globals()
        self._barcode_pair_to_int = {''.join(vv): i for i, vv in enumerate(alphabet)}
        self._int_to_barcode_pairs = {i: vv for vv, i in self._barcode_pair_to_int.items()}

    def encode(self, barcode):
        ff.set_globals(*self.ff_globals)
        barcode_as_int = [self._barcode_pair_to_int[''.join(barcode[i:i + 2])] for i in range(0, len(barcode), 2)]
        barcode_encoded_as_polynomial = self._barcode_coder.encode(barcode_as_int, return_string=False)
        barcode_encoded = ''.join([self._int_to_barcode_pairs[z] for z in barcode_encoded_as_polynomial])
        return barcode_encoded

    def decode(self, barcode_encoded):
        ff.set_globals(*self.ff_globals)
        barcode_encoded_as_int = [self._barcode_pair_to_int[''.join(barcode_encoded[i:i + 2])]
                                  for i in range(0, len(barcode_encoded), 2)]
        if self._barcode_coder.check(barcode_encoded_as_int):
            return barcode_encoded[0:self._barcode_len]
        else:
            try:
                barcode_as_int, rs_as_int = self._barcode_coder.decode(barcode_encoded_as_int, nostrip=True,
                                                                       return_string=False)
                if not self._barcode_coder.check(barcode_as_int + rs_as_int):
                    raise RSCodecError
            except RSCodecError:
                raise RSCodecError
            barcode = []
            for i in barcode_as_int:
                barcode += list(self._int_to_barcode_pairs[i])
            return barcode


# class RSPayloadAdapter:
#     def __init__(self, bits_per_z, payload_len, payload_rs_len):
#         self.bits_per_z = bits_per_z
#         self.payload_len = payload_len
#         alphabet = ['Z{}'.format(i) for i in range(1, 2**bits_per_z+1)]
#         n = payload_len + payload_rs_len
#         k = payload_len
#         c_exp = bits_per_z
#         generator = 3
#         prim = ff.find_prime_polynomials(generator=generator, c_exp=c_exp, fast_primes=False, single=True)
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


class RSPayloadAdapter:
    def __init__(self, bits_per_z, payload_len, payload_redundancy_len, payload_rs_len, vt_syndrome_n:int, vt_syndrome_k:int, k_mer_representative_to_z: Dict, z_to_k_mer_representative: Dict, binary_to_z: Dict, binary_to_k_mer_representation: Dict, bits_per_syndrome: int):
        ###########################
        self.vt_syndrome = None
        # # VT syndrome parameters for alphabet size 7
        # self.bits_per_syndrome = 3
        # payload_len = 5
        # payload_rs_len = 2
        # # VT syndrome parameters for alphabet size 16
        # bits_per_z = 4
        # payload_len = 4
        # payload_rs_len = 4
        ###########################
        self.bits_per_syndrome = bits_per_syndrome
        self.bits_per_z = bits_per_z
        self.payload_len = payload_len - payload_redundancy_len
        self.payload_rs_len = payload_rs_len
        self.payload_redundancy_len = payload_redundancy_len
        self.vt_syndrome_n = vt_syndrome_n
        self.vt_syndrome_k = vt_syndrome_k
        self.vt_syndrome = VTSyndrome(n=vt_syndrome_n, k=vt_syndrome_k,
                                      bits_per_syndrome=bits_per_syndrome,
                                      payload_len=payload_len,
                                      payload_redundancy_len=payload_redundancy_len,
                                      binary_to_k_mer_representation=binary_to_k_mer_representation)
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.binary_to_z = binary_to_z
        self.z_to_k_mer_representative = z_to_k_mer_representative
        self.binary_to_k_mer_representation = binary_to_k_mer_representation
        alphabet = ['Z{}'.format(i) for i in range(1, 2 ** bits_per_z + 1)]
        # alphabet = ['Z{}'.format(i) for i in range(1, 2 ** 6 + 1)] #TODO: change this to 2**bits_per_z with bits_per_z = 6
        n = self.payload_len + payload_rs_len
        k = self.payload_len
        # alphabet = list(range(self.vt_syndrome_n))
        c_exp = bits_per_syndrome
        generator = 2
        prim = ff.find_prime_polynomials(generator=generator, c_exp=c_exp, fast_primes=False, single=True)

        self._payload_coder_rs = rs.RSCoder(n=n, k=k, generator=generator, prim=prim, c_exp=c_exp)
        self.ff_globals = ff.get_globals()
        self._payload_to_int = {''.join(vv): i for i, vv in enumerate(alphabet)}
        self._int_to_payload = {i: vv for vv, i in self._payload_to_int.items()}

    def encode(self, payload, payload_encoded_after_vt_syndrome, binary_list, xs_array):
        ff.set_globals(*self.ff_globals)
        # syndrome_array = []
        # xs_vector_array = []
        # xs_array = []
        # binary_array = [tuple(int(char) for char in string) for string in binary_list]
        # for item in binary_array[:self.payload_len]:
        #     syn, xs_vector = self.vt_syndrome.get_syn_from_input_bits(
        #         bits_tuple=item), self.vt_syndrome.encode_message(message=item)
        #     syndrome_array.append(syn)
        #     xs_vector_array.append(xs_vector)
        #     xs_array.append(self.binary_to_k_mer_representation[item])
        converted_from_binary_to_x = []

        payload_encoded_as_polynomial = self._payload_coder_rs.encode_fast(payload_encoded_after_vt_syndrome, return_string=False)

        binary_list_info_to_concat_with_redundancy = binary_list[-self.payload_redundancy_len:]
        payload_encoded = payload[:self.payload_len]
        for binary_info in binary_list_info_to_concat_with_redundancy:
            for i,decimal_number in zip(range(0, len(binary_info), self.bits_per_syndrome),payload_encoded_as_polynomial[-self.payload_rs_len:]):
                chunk = binary_info[i:i + self.bits_per_syndrome]
                binary_array = uts.decimal_to_bits(decimal_number=decimal_number, amount_bits=self.bits_per_syndrome)
                binary_with_info_n_redundancy = uts.convert_binary_string_to_tuple(binary_array + chunk)
                xs = self.binary_to_k_mer_representation[binary_with_info_n_redundancy]
                z_vector = self.binary_to_z[binary_with_info_n_redundancy]
                payload_encoded.append(z_vector)
                xs_array.append(xs)
        return payload_encoded

    def decode(self, payload_encoded: list, erasures_pos: list):
        ff.set_globals(*self.ff_globals)
        # payload_as_int = [self._payload_to_int[z] for z in payload_encoded]
        # if self._payload_coder.check_fast(payload_as_int):
        if self._payload_coder_rs.check_fast(payload_encoded):
            # return payload_encoded[0:self.payload_len]
            return payload_encoded
        else:
            try:
                # payload_as_gf, rs_as_gf = self._payload_coder.decode(payload_as_int, nostrip=True, return_string=False)
                payload_as_gf, rs_as_gf = self._payload_coder_rs.decode(payload_encoded, erasures_pos=erasures_pos, nostrip=True, return_string=False)
            except RSCodecError:
                # return payload_encoded[0:self.payload_len]
                return payload_encoded
            # payload = [self._int_to_payload[i] for i in payload_as_gf]
            # return payload
            return payload_as_gf + rs_as_gf


class RSWideAdapter:
    def __init__(self, bits_per_z, payload_len, payload_rs_len):
        self.bits_per_z = bits_per_z
        self.payload_len = payload_len
        alphabet = ['Z{}'.format(i) for i in range(1, 2 ** bits_per_z + 1)]
        n = payload_len + payload_rs_len
        k = payload_len
        c_exp = bits_per_z
        generator = 3
        prim = ff.find_prime_polynomials(generator=generator, c_exp=c_exp, fast_primes=False, single=True)

        self._payload_coder = rs.RSCoder(n=n, k=k, generator=generator, prim=prim, c_exp=c_exp)
        self.ff_globals = ff.get_globals()
        self._payload_to_int = {''.join(vv): i for i, vv in enumerate(alphabet)}
        self._int_to_payload = {i: vv for vv, i in self._payload_to_int.items()}

    def encode(self, payload):
        ff.set_globals(*self.ff_globals)
        payload_as_int = [self._payload_to_int[z] for z in payload]
        ff.set_globals(*self.ff_globals)
        payload_encoded_as_polynomial = self._payload_coder.encode_fast(payload_as_int, return_string=False)
        payload_encoded = [self._int_to_payload[z] for z in payload_encoded_as_polynomial]
        return payload_encoded

    def decode(self, payload_encoded, erasures_pos: list) -> list:
        ff.set_globals(*self.ff_globals)
        payload_as_int = [self._payload_to_int[z] for z in payload_encoded] # TODO: origin code, lets see if the new payload_as_int code is good enough
        # payload_as_int = []
        # for z in payload_encoded:
        #     if z == 'Z0':
        #         payload_as_int.append(0)
        #     else:
        #         payload_as_int.append(self._payload_to_int[z])

        if self._payload_coder.check_fast(payload_as_int):
            return payload_encoded[0:self.payload_len]
        else:
            try:
                payload_as_gf, rs_as_gf = self._payload_coder.decode(payload_as_int, erasures_pos=erasures_pos, nostrip=True, return_string=False)
            except RSCodecError:
                return payload_encoded[0:self.payload_len]
            payload = [self._int_to_payload[i] for i in payload_as_gf]
            return payload
