import itertools
from typing import Dict

import numpy as np
from unireedsolomon.unireedsolomon import rs, RSCodecError
from unireedsolomon.unireedsolomon import ff
from VTSyndrom import VTSyndrom


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
    def __init__(self, bits_per_z, payload_len, payload_rs_len, vt_syndrome_n:int, vt_syndrome_k:int, k_mer_representative_to_z: Dict):
        ###########################
        self.vt_syndrome = None
        bits_per_z = 3
        payload_len = 4
        payload_rs_len = 3
        ###########################
        self.bits_per_z = bits_per_z
        self.payload_len = payload_len
        self.vt_syndrome_n = vt_syndrome_n
        self.vt_syndrome_k = vt_syndrome_k
        self.init_vt_syndrome(n=self.vt_syndrome_n, k=self.vt_syndrome_k)
        self.k_mer_representative_to_z = k_mer_representative_to_z
        # alphabet = ['Z{}'.format(i) for i in range(1, 2 ** bits_per_z + 1)]
        alphabet = ['Z{}'.format(i) for i in range(1, 2 ** 6 + 1)] #TODO: change this to 2**bits_per_z with bits_per_z = 6
        n = payload_len + payload_rs_len
        k = payload_len
        # alphabet = list(range(self.vt_syndrome_n))
        c_exp = bits_per_z
        generator = 2
        prim = ff.find_prime_polynomials(generator=generator, c_exp=c_exp, fast_primes=False, single=True)

        self._payload_coder = rs.RSCoder(n=n, k=k, generator=generator, prim=prim, c_exp=c_exp)
        self.ff_globals = ff.get_globals()
        self._payload_to_int = {''.join(vv): i for i, vv in enumerate(alphabet)}
        self._int_to_payload = {i: vv for vv, i in self._payload_to_int.items()}

    def encode(self, payload, binary_list):
        ff.set_globals(*self.ff_globals)
        syndrome_array = []
        xs_vector_array = []
        binary_array = [[int(char) for char in string] for string in binary_list]
        for item in binary_array:
            syn, ind, xs_vector = self.vt_syndrome.one_letter(self.table, item)
            syndrome_array.append(syn)
            xs_vector_array.append(xs_vector)
        converted_from_binary_to_x = []

        payload_encoded_as_polynomial = self._payload_coder.encode_fast(syndrome_array, return_string=False)

        for binary_array in xs_vector_array:
            xi_values = [f'X{i + 1}' for i, bit in enumerate(binary_array) if bit == 1]
            xi_values = tuple(xi_values)
            values = [f'{i + 1}' for i, bit in enumerate(binary_array) if bit == 1]
            converted_from_binary_to_x.append(xi_values)
        # payload_as_int = [self._payload_to_int[z] for z in payload]
        # payload_encoded_as_polynomial = self._payload_coder.encode_fast(payload_as_int, return_string=False)
        z_array = [self._payload_to_int[self.k_mer_representative_to_z[values]] for values in converted_from_binary_to_x]
        payload_encoded_as_polynomial = self._payload_coder.encode_fast(z_array, return_string=False)
        payload_encoded = [f'X{i}' for i in payload_encoded_as_polynomial]
        return payload_encoded

    def decode(self, payload_encoded):
        ff.set_globals(*self.ff_globals)
        # payload_as_int = [self._payload_to_int[z] for z in payload_encoded]
        # if self._payload_coder.check_fast(payload_as_int):
        if self._payload_coder.check_fast(payload_encoded):
            return payload_encoded[0:self.payload_len]
        else:
            try:
                # payload_as_gf, rs_as_gf = self._payload_coder.decode(payload_as_int, nostrip=True, return_string=False)
                payload_as_gf, rs_as_gf = self._payload_coder.decode(payload, nostrip=True, return_string=False)
            except RSCodecError:
                return payload_encoded[0:self.payload_len]
            payload = [self._int_to_payload[i] for i in payload_as_gf]
            return payload

    def init_vt_syndrome(self, n, k):
        self.vt_syndrome = VTSyndrom(n=n, k=k)
        self.table = self.vt_syndrome.create_table()


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
