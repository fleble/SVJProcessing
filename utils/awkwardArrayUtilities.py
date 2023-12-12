import awkward as ak


def as_type(ak_array, dtype):
    if ak.count(ak_array) == 0:
        dummy_ak_array = ak.values_astype(ak.Array([[0.]]), dtype)
        ak_array = ak.concatenate((ak_array, dummy_ak_array), axis=0)
        ak_array = ak_array[:-1]
    else:
        ak_array = ak.values_astype(ak.from_iter(ak_array), dtype)

    return ak_array

