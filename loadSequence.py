from Bio import SeqIO
# For more info: https://biopython.org/wiki/SeqIO
from operator import attrgetter
import pickle
from models import Kmer, Hit, FragmentExtension

ALP = "ACGTN"


def load_seq(sequence_file_name):
    print('\n\n~ load_seq ~ fasta_file: ' + sequence_file_name + '\n')
    # The argument 'rU' means open for reading using universal readline mode
    # To handle a list of records:
    # records = list(SeqIO.parse("example.fasta", "fasta"))

    # For Python 2.7 use "rU"
    with open(sequence_file_name, "r") as sequence_file:
        for record in SeqIO.parse(sequence_file, "fasta"):
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
    print('\n\n~ compute_frecuency ~ sequence: ' + sequence + '\n')
    frequency = [0] * 4 ** kmer_size
    i = 0
    while i < len(sequence) - (kmer_size - 1):
        n = kmer_index(sequence[i:i + kmer_size], kmer_size, alphabet)
        frequency[n] += 1
        i += 1
        # print("word: %s, freq: %d" % (kmer_index2word(n, 5, alphabet), frequency[n]))
    return frequency


def generate_dictionary(sequence, kmer_size, alphabet, o_file_name):
    print('\n\n~ generate_dictionary ~ output: ' + o_file_name + '\n')
    kmers = []
    position = 0

    # Index sequence and generate dictionary
    while position < len(sequence) - (kmer_size - 1):
        kmer_hash = kmer_index(sequence[position:position + kmer_size], kmer_size, alphabet)
        kmer = Kmer(kmer_hash, position)
        kmers.append(kmer)
        position += 1
    # Sort by kmer hash
    sorted_kmers = sorted(kmers, key=attrgetter('k_hash'))
    # sorted_kmers = sorted(kmers, key=lambda k: k.k_hash)

    # Group by kamer hash
    # for i in sorted_kmers:
    #   print(str(i))
    grouped_kmers = group_dictionary(sorted_kmers)

    # Write sorted array as plain text to file
    #   with open(file_name, 'w') as seq_dic:
    #       for item in sorted_kmers:
    #           seq_dic.write('%s\n' % item)

    # Dump array into file using pickle
    with open(o_file_name, 'wb') as seq_dic:
        pickle.dump(grouped_kmers, seq_dic)
    seq_dic.close()


def group_dictionary(dic):
    print('\n\n~ group_dictionary ~\n')
    grouped_dic = []
    positions = []
    i = 1
    current_kmer = dic[0]
    prep_kmer = current_kmer
    positions.append(current_kmer.position)
    while i < len(dic):
        if current_kmer.k_hash == dic[i].k_hash:
            positions.append(dic[i].position)
            i += 1
        else:
            k = Kmer(current_kmer.k_hash, positions)
            grouped_dic.append(k)
            #print(kmer_index2word(k.k_hash, 5, ALP), str(k))
            positions = []

            prep_kmer = current_kmer
            current_kmer = dic[i]
            positions.append(current_kmer.position)
            i += 1
    # Check last kmer
    if prep_kmer.k_hash != current_kmer.k_hash:
        k = Kmer(current_kmer.k_hash, positions)
        grouped_dic.append(k)
        # print(kmer_index2word(k.k_hash, 5, ALP), str(k))
    return grouped_dic


def generate_hits(i_file_name_1, i_file_name_2, o_file_name):
    print('\n\n~ generate_hits ~ dic_1: ' + i_file_name_1 + ', dic_2: ' + i_file_name_2 + ', output: ' + o_file_name + '\n')
    hits = []
    with open(i_file_name_1, 'rb') as dic_file_1, \
            open(i_file_name_2, 'rb') as dic_file_2, \
            open(o_file_name, 'wb') as hits_file:

        # Load both sorted dictionaries
        dic_1 = pickle.load(dic_file_1)
        dic_2 = pickle.load(dic_file_2)
        i = 0
        j = 0
        while i < len(dic_1) and j < len(dic_2):
            if dic_1[i].k_hash == dic_2[j].k_hash:
                for p_1 in dic_1[i].position:
                    for p_2 in dic_2[j].position:
                        hit = Hit(dic_1[i].k_hash, p_1, p_2)
                        hits.append(hit)
                        #print(str(hit))
                        #print(kmer_index2word(hit.k_hash, 5, ALP), str(hit))
                # print(str(hit))
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


# alineamiento mediante fragmentos sin huecos, score por longitud de secuencia
def fragment_extension_from_hits(hits_file_name, output_file_name, sequence_1, sequence_2, kmer_size):
    with open(hits_file_name, 'rb') as hits_file, open(output_file_name, 'wb') as output_file:
        hits = pickle.load(hits_file)
        results = []
        for hit in hits:
            # Tener en cuenta la matriz PAM
            left_score = kmer_size
            right_score = kmer_size

            max_left_score = kmer_size
            max_right_score = kmer_size
            num_identity = kmer_size


            # contamos por la izquierda
            left_offset = 1
            right_offset = 1

            print 'mirando por la izquierda'
            while left_score >= 0 and \
                    (hit.position_1 - left_offset) >= 0 and \
                    (hit.position_2 - left_offset) >= 0:

                if sequence_1[hit.position_1 - left_offset] == sequence_2[hit.position_2 - left_offset]:
                    left_score += 1
                    num_identity += 1
                    if left_score > max_left_score:
                        max_left_score = left_score
                else:
                    left_score -= 1

                #'''
                print('S1.x_start',
                      hit.position_1,
                      hit.position_1 - left_offset,
                      sequence_1[hit.position_1 - left_offset],
                      left_score,
                      max_left_score,
                      num_identity)
                print('S2.y_start',
                      hit.position_2,
                      hit.position_2 - left_offset,
                      sequence_2[hit.position_2 - left_offset],
                      left_score,
                      max_left_score,
                      num_identity)
                #'''
                left_offset += 1

            # comprobar que no se salga de los limites
            if (hit.position_1 - left_offset) >= 0 and (hit.position_2 - left_offset) >= 0:
                x_start = hit.position_1 - left_offset  # fragmento mas largo, sin tener en cuenta el score
                y_start = hit.position_2 - left_offset
            else:
                x_start = hit.position_1  # fragmento mas largo, sin tener en cuenta el score
                y_start = hit.position_2
            # print(x_start, y_start)

            # contamos por la derecha
            print 'mirando por la derecha'
            while right_score >= 0 and \
                    (hit.position_1 + kmer_size + right_offset) < len(sequence_1) and \
                    (hit.position_2 + kmer_size + right_offset) < len(sequence_2):

                if sequence_1[hit.position_1 + kmer_size + right_offset] == sequence_2[hit.position_2 + kmer_size + right_offset]:
                    right_score += 1
                    num_identity += 1
                    if right_score > max_right_score:
                        max_right_score = right_score
                else:
                    right_score -= 1

                #'''
                print('S1.',
                      hit.position_1 + kmer_size,
                      hit.position_1 + kmer_size + right_offset,
                      sequence_1[hit.position_1 + kmer_size + right_offset],
                      right_score,
                      max_right_score,
                      num_identity)
                print('S2.',
                      hit.position_2 + kmer_size,
                      hit.position_2 + kmer_size + right_offset,
                      sequence_2[hit.position_2 + kmer_size + right_offset],
                      right_score,
                      max_right_score,
                      num_identity)
                #'''
                right_offset += 1

            # comprobar que no se salga de los limites
            if (hit.position_1 + kmer_size + right_offset) < len(sequence_1) and \
                    (hit.position_2 + kmer_size + right_offset) < len(sequence_2):
                x_end = hit.position_1 + kmer_size + right_offset
                y_end = hit.position_2 + kmer_size + right_offset
            else:
                x_end = hit.position_1 + kmer_size
                y_end = hit.position_2 + kmer_size
                # print(x_end, y_end)

            # comprobar que no se haya metido en el bucle
            # print(x_start, x_end, y_start, y_end)
            if x_start >= 0 and y_start >= 0 and x_end <= len(sequence_1) and y_end <= len(sequence_1):
                score = max_left_score + max_right_score - kmer_size
                fragment = FragmentExtension(x_start, x_end, y_start, y_end, score, num_identity)
                results.append(fragment)
                print(str(fragment))

        pickle.dump(results, output_file)


# alineamiento mediante fragmentos sin huecos, mediante score mayor
def fragment_extension_from_hits2(hits_file_name, output_file_name, sequence_1, sequence_2, kmer_size):
    print('\n\n~ fragment_extension_from_hits2 ~ hits_file: ' + hits_file_name + ', output: ' + output_file_name + '\n')
    with open(hits_file_name, 'rb') as hits_file, open(output_file_name, 'wb') as output_file:
        hits = pickle.load(hits_file)
        results = []
        print(len(hits))
        i = 0
        for hit in hits:
            # Tener en cuenta la matriz PAM
            left_score = kmer_size
            right_score = kmer_size

            max_left_score = kmer_size
            max_right_score = kmer_size

            best_x_start = hit.position_1
            best_y_start = hit.position_2

            best_x_end = hit.position_1 + kmer_size
            best_y_end = hit.position_1 + kmer_size

            num_identity = kmer_size

            # contamos por la izquierda
            left_offset = 1
            right_offset = 1

            # print 'mirando por la izquierda'
            while left_score >= 0 and \
                    (hit.position_1 - left_offset) >= 0 and \
                    (hit.position_2 - left_offset) >= 0:

                if sequence_1[hit.position_1 - left_offset] == sequence_2[hit.position_2 - left_offset]:
                    left_score += 1
                    num_identity += 1
                    if left_score >= max_left_score:
                        max_left_score = left_score
                        best_x_start = hit.position_1 - left_offset
                        best_y_start = hit.position_2 - left_offset
                        # print(best_x_start, best_y_start)
                else:
                    left_score -= 1
                '''
                print '==='
                print('S1.x_start',
                      hit.position_1,
                      hit.position_1 - left_offset,
                      sequence_1[hit.position_1 - left_offset],
                      left_score,
                      max_left_score,
                      num_identity)
                print(best_x_start, best_y_start)
                print('S2.y_start',
                      hit.position_2,
                      hit.position_2 - left_offset,
                      sequence_2[hit.position_2 - left_offset],
                      left_score,
                      max_left_score,
                      num_identity)
                '''
                left_offset += 1

            # contamos por la derecha
            # print 'mirando por la derecha'
            while right_score >= 0 and \
                    (hit.position_1 + kmer_size + right_offset) < len(sequence_1) and \
                    (hit.position_2 + kmer_size + right_offset) < len(sequence_2):

                if sequence_1[hit.position_1 + kmer_size + right_offset] == sequence_2[hit.position_2 + kmer_size + right_offset]:
                    right_score += 1
                    num_identity += 1
                    if right_score >= max_right_score:
                        max_right_score = right_score
                        best_x_end = hit.position_1 + kmer_size + right_offset
                        best_y_end = hit.position_2 + kmer_size + right_offset
                else:
                    right_score -= 1
                '''
                print '==='
                print('S1.',
                      hit.position_1 + kmer_size,
                      hit.position_1 + kmer_size + right_offset,
                      sequence_1[hit.position_1 + kmer_size + right_offset],
                      right_score,
                      max_right_score,
                      num_identity)
                print(best_x_end, best_y_end)
                print('S2.',
                      hit.position_2 + kmer_size,
                      hit.position_2 + kmer_size + right_offset,
                      sequence_2[hit.position_2 + kmer_size + right_offset],
                      right_score,
                      max_right_score,
                      num_identity)
                '''
                right_offset += 1

            # comprobar que no se haya metido en el bucle
            # print(x_start, x_end, y_start, y_end)
            # if x_start >= 0 and y_start >= 0 and x_end <= len(sequence_1) and y_end <= len(sequence_1):
            score = max_left_score + max_right_score - kmer_size
            fragment = FragmentExtension(best_x_start, best_x_end, best_y_start, best_y_end, score, num_identity)
            results.append(fragment)
            # print(str(fragment))
            i += 1

            if i % 1000 == 0:
                print('Llevo ' + str(i) + ' de ' + str(len(hits)))

        pickle.dump(results, output_file)


# alineamiento mediante fragmentos sin huecos, mediante score mayor
def fragment_extension_from_hits_to_csv(hits_file_name, output_file_name, sequence_1, sequence_2, kmer_size):
    print('\n\n~ fragment_extension_from_hits2 ~ hits_file: ' + hits_file_name + ', output: ' + output_file_name + '\n')
    with open(hits_file_name, 'rb') as hits_file, open(output_file_name, 'w') as output_file:
        print('Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY\n')
        hits = pickle.load(hits_file)
        # results = []
        print('Total hits: ' + str(len(hits)))
        i = 0
        for hit in hits:
            # Tener en cuenta la matriz PAM
            left_score = kmer_size
            right_score = kmer_size

            max_left_score = kmer_size
            max_right_score = kmer_size

            best_x_start = hit.position_1
            best_y_start = hit.position_2

            best_x_end = hit.position_1 + kmer_size
            best_y_end = hit.position_2 + kmer_size

            num_identity = kmer_size
            max_identity = kmer_size

            # contamos por la izquierda
            left_offset = 1
            right_offset = 1

            #print 'mirando por la izquierda'
            while left_score >= 0 and \
                    (hit.position_1 - left_offset) >= 0 and \
                    (hit.position_2 - left_offset) >= 0:

                if sequence_1[hit.position_1 - left_offset] == sequence_2[hit.position_2 - left_offset]:
                    left_score += 1
                    num_identity += 1
                    if left_score >= max_left_score:
                        max_left_score = left_score
                        best_x_start = hit.position_1 - left_offset
                        best_y_start = hit.position_2 - left_offset
                        # print(best_x_start, best_y_start)
                        max_identity = num_identity
                else:
                    left_score -= 1
                left_offset += 1

            # contamos por la derecha
            #print 'mirando por la derecha'
            while right_score >= 0 and \
                    (hit.position_1 + kmer_size + right_offset) < len(sequence_1) and \
                    (hit.position_2 + kmer_size + right_offset) < len(sequence_2):

                if sequence_1[hit.position_1 + kmer_size + right_offset] == sequence_2[hit.position_2 + kmer_size +
                                                                                       right_offset]:
                    right_score += 1
                    num_identity += 1
                    if right_score >= max_right_score:
                        max_right_score = right_score
                        best_x_end = hit.position_1 + kmer_size + right_offset
                        best_y_end = hit.position_2 + kmer_size + right_offset
                        max_identity = num_identity
                else:
                    right_score -= 1
                right_offset += 1

            # Calculo de la semejanza
            score = max_left_score + max_right_score - kmer_size
            length_fragment = best_x_end - best_x_start
            # length = best_y_end - best_y_start
            similarity = round(float(score) / length_fragment,2)
            identity_percentage = round(float(max_identity * 100) / length_fragment,2)

            # Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY
            # if identity_percentage > 50:
            '''output_file.write('Frag,' + str(best_x_start) + ',' + str(best_y_start) + ',' + str(best_x_end) + ',' +
                              str(best_y_end) + ',f,0,' + str(length_fragment) + ',' + str(score) + ',' +
                              str(max_identity) + ',' + str(similarity) + ',' + str(identity_percentage) + ',0,0\n')
            '''
            output_file.write('Frag,' + str(best_x_start) + ',' + str(best_y_start) + ',' + str(best_x_end) + ',' +
                              str(best_y_end) + ',f,0,' + str(length_fragment) + ',' + str(score) + ',' +
                              str(max_identity) + ',' + str(similarity) + ',' + str(identity_percentage) + ',0,0\n')
            # fragment = FragmentExtension(best_x_start, best_x_end, best_y_start, best_y_end, score, num_identity)
            # results.append(fragment)
            # print(str(fragment))

            '''print('final left: ' + str(left_score) + ', pos_izq:' + str(hit.position_1 - left_offset) + ', final right: '
                 + str(right_score) + ', pos_dcha: ' + str(hit.position_1 + kmer_size + right_offset))
            
            print '======='
            '''
            i += 1
            #print score, length_fragment, similarity

            if i % 1000 == 0:
                print(str(i) + ' computed')

        #pickle.dump(results, output_file)


# -------- WORKFLOW --------

#K_SIZE = 5
# s1 = load_seq("example.fasta")
# s2 = load_seq("example2.fasta")

# For Python 3 print is a function
# print(kmer_index(s1, K_SIZE, ALP))
# print(kmer_index2word(75, K_SIZE, ALP))
# print(compute_frequency(s1, K_SIZE, ALP))

#generate_dictionary(s1, K_SIZE, ALP, 'seq_dic_1')
#generate_dictionary(s2, K_SIZE, ALP, 'seq_dic_2')
#generate_hits('seq_dic_1', 'seq_dic_2', 'o_hits')
#fragment_extension_from_hits2('o_hits', 'o_fragment_extension', s1, s2, K_SIZE)

K_SIZE = 20

s1 = load_seq("mycoplasma_hyoneumoniae_232_simple.fasta")
s2 = load_seq("mycoplasma_hyoneumoniae_7422_simple.fasta")

# generate_dictionary(s1, K_SIZE, ALP, '232_dic_1')
# generate_dictionary(s2, K_SIZE, ALP, '7422_dic_2')
# generate_hits('232_dic_1', '7422_dic_2', 'mh_hits')
# fragment_extension_from_hits2('mh_hits', 'mh_fragment_extension', s1, s2, K_SIZE)
fragment_extension_from_hits_to_csv('mh_hits', 'mh_fragment_extension.csv', s1, s2, K_SIZE)

