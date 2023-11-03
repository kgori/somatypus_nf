class SomatypusError(Exception):
    pass


class Interval():
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper

    def __contains__(self, item):
        return self.lower <= item <= self.upper

    def __repr__(self):
        return 'Interval({}, {})'.format(self.lower, self.upper)

    def __lt__(self, other):
        return (self.lower, self.upper) < (other.lower, other.upper)

    def __eq__(self, other):
        return (self.lower, self.upper) == (other.lower, other.upper)

    def intersects(self, other):
        return other.lower <= self.upper and other.upper >= self.lower

    def merge(self, other):
        if not self.intersects(other):
            raise ValueError('Intervals do not intersect')
        return Interval(min(self.lower, other.lower), max(self.upper, other.upper))

    def absorb(self, other):
        if self.intersects(other):
            self.lower = min(self.lower, other.lower)
            self.upper = max(self.upper, other.upper)
            return True
        return False


class Variant():
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        try:
            self.pos = int(pos)
        except ValueError:
            raise SomatypusError('Position is not an integer: {}'.format(pos))
        self.ref = ref
        self.alt = alt

    def __repr__(self):
        return '{}:{},{}>{}'.format(self.chrom, self.pos, self.ref, self.alt)

    def __lt__(self, other):
        return (self.chrom, self.pos, self.ref, self.alt) < (other.chrom, other.pos, other.ref, other.alt)

    def __eq__(self, other):
        return (self.chrom, self.pos, self.ref, self.alt) == (other.chrom, other.pos, other.ref, other.alt)

    def __hash__(self):
        return hash((self.chrom, self.pos, self.ref, self.alt))

    def is_snv(self):
        return len(self.ref) == len(self.alt) == 1

    def is_del(self):
        return len(self.ref) > len(self.alt)

    def is_ins(self):
        return len(self.ref) < len(self.alt)

    def is_indel(self):
        return self.is_del() or self.is_ins()

    def is_multiallelic(self):
        return ',' in self.alt


def get_variant_from_line(line):
    col = line.strip().split('\t')
    chrom = col[0]
    pos = col[1]
    ref = col[3]
    alt = col[4]
    return Variant(chrom, pos, ref, alt)
