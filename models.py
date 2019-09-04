
class Kmer:
    def __init__(self, k_hash, position):
        self.k_hash = k_hash
        self.position = position

    def __str__(self):
        return u'%s,%s' % (str(self.k_hash), str(self.position))


class Hit:
    def __init__(self, k_hash, position_1, position_2):
        self.k_hash = k_hash
        self.position_1 = position_1
        self.position_2 = position_2

    def __str__(self):
        return u'%s,%s,%s' % (str(self.k_hash), str(self.position_1), str(self.position_2))


class FragmentExtension:
    def __init__(self, x_start, x_end, y_start, y_end, score, num_indentity):
        self.x_start = x_start
        self.x_end = x_end
        self.y_start = y_start
        self.y_end = y_end
        self.score = score
        self.num_indentity = num_indentity

    def __str__(self):
        return u'%s,%s | %s,%s | %s,%s' % (str(self.x_start), str(self.x_end), str(self.y_start), str(self.y_end),
                                       str(self.score), str(self.num_indentity))
