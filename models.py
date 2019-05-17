
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
