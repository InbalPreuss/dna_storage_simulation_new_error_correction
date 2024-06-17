from itertools import combinations
import numpy as np
class VTSyndrom:
    def __init__(self, n, k):
        self.n = n
        self.k = k
    def bits_to_number(self, bits):
        number = 0
        for i, bit in enumerate(reversed(bits)):  # reverse to start from LSB
            number += bit * (2 ** i)
        return number
    def create_word(self, ones_positions):
        bits = [0] * self.n
        for pos in ones_positions:
            bits[pos] = 1
        return bits

    def generate_table(self):
        table = [[],[],[],[],[],[],[],[]]
        for ones_positions in combinations(range(self.n), self.k):
            bits = self.create_word(ones_positions, self.n)
            syndrome_sum = sum((i + 1) * bit for i, bit in enumerate(bits))
            #sum = np.sum((1 + np.arange(n_y)) * y)
            syndrom_index = np.mod(syndrome_sum, self.n)
            if len(table[syndrom_index]) < 8:
                table[syndrom_index].append(bits)
        return table
    def print_table(self, table):
        for i in range(8):
            print(table[i])
            print(len(table[i]))


    def one_letter(self, table, my_input):
         #input 6 bits
        #Step 1 - convert input to sigma
        print("input is ", self.bits_to_number(my_input))
        #Step  2 -
        syndrom = self.bits_to_number(my_input[0:3])
        index = self.bits_to_number(my_input[3:])
        print("syndrom = ", syndrom)
        print("index = ", index)
        return syndrom, index, table[syndrom][index]

    def create_table(self):
         table = self.generate_table(self.n, self.k)
         return table
         # chosen_vector = table[syndrom, index]#4 zeros and 4 ones
         # print_table(table)

"""
def main():
    n = 8
    k = 4
    vtsynd = VTSyndrom(n,k)
    table = vtsynd.create_table()
    #xs_vector is the 8 bits where 4 ones and 4 zeros
    #run the code on 6*6 input to create 6 letters
    syn, ind, xs_vector =vtsynd.one_letter(table, [0,0,1,0,0,1])
    print(syn)
    vtsynd.one_letter(table, [0, 0, 1, 0, 0, 1])
    vtsynd.one_letter(table, [0, 0, 1, 0, 0, 1])
    vtsynd.one_letter(table, [0, 0, 1, 0, 0, 1])
    vtsynd.one_letter(table, [0, 0, 1, 0, 0, 1])

main()
"""