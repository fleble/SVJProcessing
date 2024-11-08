import awkward as ak


def as_type(ak_array, dtype):
    if ak.count(ak_array) == 0:
        dummy_ak_array = ak.values_astype(ak.Array([[0.]]), dtype)
        ak_array = ak.concatenate((ak_array, dummy_ak_array), axis=0)
        ak_array = ak_array[:-1]
    else:
        ak_array = ak.values_astype(ak.from_iter(ak_array), dtype)

    return ak_array


        
def divide_ak_arrays(ak_array1, ak_array2, division_by_zero_value=1., verbose=False):
    """Makes the division of an ak array by another one.
    
    The arrays ak_array1 and ak_array2 must have the same jagged structure.
    If division by zero for some indices, a default value to use can be
    defined, see examples.
    The output array has the same jagged structure as the input ak arrays.
    

    Args:
        ak_array1 (awkward.Array[float])
        ak_array2 (awkward.Array[float]) 
        division_by_zero_value (float, optional, default=1.)
        verbose (bool, optional, default=False)

    Returns:
        awkward.Array[float]: ak_array1 / ak_array2

    Examples:
        >>> ak_array1 = ak.Array([ [0, 3], [5], [1] ])
        >>> ak_array2 = ak.Array([ [3, 3], [0], [2] ])
        >>> divide_ak_arrays(ak_array1, ak_array2)
        [ [0, 1], [1], [0.5] ]
    """

    is_not_zero = (ak_array2!=0.)
    if (not ak.all(is_not_zero)) and verbose:
        print("The following warning about true_divide can be safely ignored.")

    raw_division = ak_array1/ak_array2
    division = ak.where(is_not_zero, raw_division, division_by_zero_value*ak.ones_like(ak_array1))

    # This implementation seems slower:
    #division = ak.Array([ [ x1/x2 if x2 != 0. else division_by_zero_value for x1, x2 in zip(y1, y2) ] for y1, y2 in zip(ak_array1, ak_a0rray2) ])

    return division



def is_in(array1, array2):
    """
    Args:
        array1 (awkward.highlevel.Array[int])
        array2 (awkward.highlevel.Array[int])

    Returns:
        awkward.highlevel.Array[bool]

    Examples:
        >>> example_array1 = [[0, 1, 1, 2, 2, 3], [0, 0, 1], [0, 1, 1, 2]]
        >>> example_array2 = [[0, 3], [], [0, 2]]
        >>> is_in(example_array1, example_array2)
        [[True, False, False, False, False, True], [False, False, False] [True, False, False, True]]
    """

    return ak.Array([[True if x in y2 else False for x in y1] for y1, y2 in zip(array1, array2)])


def is_in_list(ak_array, list_):
    """Check whether the elements of an ak array are in a list of elements.
    
    Args:
        ak_array1 (awkward.Array[T])
        list_ (list[T])

    Returns:
        awkward.Array[bool]

    Examples:
        >>> ak_array = ak.Array([ [11, 22, -11], [22], [111, 211, -11, 11] ])
        >>> list_ = [11, -11]
        >>> is_in(ak_array, list_)
        [[True, False, True], [False], [False, False, True, True]]
    """

    ak_bool = False * ak.ones_like(ak_array, dtype=bool)
    for el in list_:
        ak_bool = ak_bool | (ak_array == el)

    return ak_bool


def sort_array_with_fields(ak_array, field_name, ascending=False):
    """Sort ak array of records using the values in one of its fields.

    Args:
        ak_array (ak.Array): Ak array of records
        field_name (str): the field to use to sort the array
        ascending (bool, optional, default=False): Set to False to sort objects
            by pT like usually done in HEP (highest pT first)
    """

    sorted_indices = ak.argsort(ak_array[field_name], ascending=ascending)
    return ak_array[sorted_indices]

def get_type(ak_array):
    return str(ak.type(ak_array)).split("*")[-1].replace("?", "").replace(" ", "")