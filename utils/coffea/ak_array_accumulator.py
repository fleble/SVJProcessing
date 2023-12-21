import awkward as ak
from coffea import processor


class AkArrayAccumulator(processor.AccumulatorABC):
    def __init__(self, value=ak.Array({})):
        self.value = value

    def identity(self):
        return AkArrayAccumulator(ak.Array({}))
    
    def add(self, other):
        if len(other.value) == 0:
            return
        elif len(self.value) == 0:
            self.value = other.value
        else:
            self.value = ak.concatenate((self.value, other.value), axis=0)

