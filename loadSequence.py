from operator import attrgetter
from Bio import SeqIO
import pickle
from models import Kmer, Hit
# For more info: https://biopython.org/wiki/SeqIO

ALP = "ACGT"


def load_seq(sequence_file):
    # The argument 'rU' means open for reading using universal readline mode
    # To handle a list of records:
    # records = list(SeqIO.parse("example.fasta", "fasta"))

    # For Python 2.7 use "rU"
    with open(sequence_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return record.seq


def change_char(string, position, replacement_char):
    return string[:position] + replacement_char + string[position + 1:]


def kmer_index(sequence, kmer_size, alphabet):
    res = 0
    if len(sequence) < kmer_size:
        return -1
    for i in range(0, kmer_size):
        # print("%s * 4^%d-%d-1: %d" % (sequence[i], kmer_size, i, 4 ** (kmer_size - i - 1)))
        res += alphabet.index(sequence[i]) * 4 ** (kmer_size - i - 1)
    return res


def kmer_index2word(kmer_hash, kmer_size, alphabet):
    quotient = kmer_hash
    word = 'A' * kmer_size
    i = 0

    while quotient >= 4:
        reminder = quotient % 4  # get letter position on the alphabet
        quotient = quotient / 4  # remove power of 4

        # For Python 3 we need to convert to int the reminder (it's calculated as a float)
        rem_int = int(reminder)

        word = change_char(word, kmer_size - i - 1, alphabet[rem_int])
        i += 1

    # For Python 3 we need to convert to int the reminder (it's calculated as a float)
    quo_int = int(quotient)
    word = change_char(word, kmer_size - i - 1, alphabet[quo_int])
    return word


def compute_frequency(sequence, kmer_size, alphabet):
    frequency = [0] * 4 ** kmer_size
    i = 0
    while i < len(sequence) - (kmer_size - 1):
        n = kmer_index(sequence[i:i + kmer_size], kmer_size, alphabet)
        frequency[n] += 1
        i += 1
        # print("word: %s, freq: %d" % (kmer_index2word(n, 5, alphabet), frequency[n]))
    return frequency


def generate_dictionary(sequence, kmer_size, alphabet, o_file_name):
    kmers = []
    position = 0
    # Index sequence and generate dictionary
    while position < len(sequence) - (kmer_size - 1):
        kmer_hash = kmer_index(sequence[position:position + kmer_size], kmer_size, alphabet)
        position += 1
        kmer = Kmer(kmer_hash, position)
        kmers.append(kmer)
    # Sort by kmer hash
    
    # Dump array into file using pickle
    with open(o_file_name, 'wb') as seq_dic:
        pickle.dump(sorted_kmers, seq_dic)
    seq_dic.close()


def generate_hits(i_file_name_1, i_file_name_2, o_file_name_3):
    hits = []
    with open(i_file_name_1, 'rb') as dic_file_1, \
            open(i_file_name_2, 'rb') as dic_file_2, \
            open(o_file_name_3, 'wb') as hits_file:
        # Load both sorted dictionaries
        dic_1 = pickle.load(dic_file_1)
        dic_2 = pickle.load(dic_file_2)
        i = 0
        j = 0
        while i < len(dic_1) and j < len(dic_2):
            if dic_1[i].k_hash == dic_2[j].k_hash:
                hit = Hit(dic_1[i].k_hash, dic_1[i].position, dic_2[j].position)
                hits.append(hit)
                print(str(hit))
                i += 1
            elif dic_1[i].k_hash < dic_2[j].k_hash:
                i += 1
            elif dic_1[i].k_hash > dic_2[j].k_hash:
                j += 1
        # Dump array into file using pickle
        pickle.dump(hits, hits_file)
    dic_file_1.close()
    dic_file_2.close()
    hits_file.close()


s1 = load_seq("example.fasta")
s2 = load_seq("example2.fasta")

# For Python 2.7
# print compute_frequency(s, 5, ALP)
# print kmer_index(s, 5, ALP)
# print kmer_index2word(75, 5, ALP)

# For Python 3 print is a function
# print(kmer_index(s1, 5, ALP))
# print(kmer_index2word(75, 5, ALP))
# print(compute_frequency(s1, 5, ALP))

generate_dictionary(s1, 5, ALP, 'seq_dic_1')
generate_dictionary(s2, 5, ALP, 'seq_dic_2')
generate_hits('seq_dic_1', 'seq_dic_2', 'o_hits')
