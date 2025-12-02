from coffea import processor


class DictAccumulator(processor.AccumulatorABC):
    def __init__(self, value={}):
        self.value = value

    def identity(self):
        return DictAccumulator({})
    
    def add(self, other):
        if len(other.value.keys()) == 0:
            return
        elif len(self.value.keys()) == 0:
            self.value = other.value
        else:
            for key in other.value.keys():
                if key in self.value.keys(): self.value[key] += other.value[key]
                else: self.value[key] = other.value[key]

