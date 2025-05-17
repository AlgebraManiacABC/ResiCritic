

class Ratio:
    '''
    Defines a pair of numbers as numerator and denominator of a ratio.
    .numer: numerator
    .denom: denominator
    .ratio(): returns numer / denom as a float
    '''
    def __init__(self, numer: int, denom: int):
        self.numer = numer
        self.denom = denom

    def ratio(self) -> float:
        return float(self.numer / self.denom) if self.denom else 'nan'

    def __str__(self):
        return f"{self.numer}:{self.denom} ({self.ratio()})"
