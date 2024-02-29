import numpy as np
import awkward as ak
import uproot

from utils import awkward_array_utilities as akUtl
from utils.Logger import *


def __cast_unknown_type_branches(ak_array, field):

    type_text = str(ak.type(ak_array))

    if "unknown" in type_text:
        log.warning(f"Branch {field} has unknown type, converting to float64")
        ak_array = akUtl.as_type(ak_array, np.float64)
        
    elif "?" in type_text:
        new_type_text = type_text.split("*")[-1].replace("?", "").replace(" ", "")
        log.warning(f"Branch {field} has unclear type {type_text}, converting to {new_type_text}")
        ak_array = akUtl.as_type(ak_array, getattr(np, new_type_text))
        
    return ak_array


def __prepare_array(ak_array, field):
    return  ak.packed(ak.without_parameters(__cast_unknown_type_branches(ak_array, field)))


def __is_rootcompat(a):
    """Is it a flat or 1-d jagged array?"""

    t = ak.type(a)
    if isinstance(t, ak._ext.ArrayType):
        if isinstance(t.type, ak._ext.PrimitiveType):
            return True
        if isinstance(t.type, ak._ext.ListType) and isinstance(
            t.type.type, ak._ext.PrimitiveType
        ):
            return True
    return False


def __zip_composite(ak_array):
    """Zip together branches belonging to the same collection.
    
    Args:
        ak_array (ak.Array): 
    
    Returns:
        ak.Array
    """

    # Additional naming scheme to allow composite object readback
    _rename_lookup = {
        "pt": "/.fPt",
        "eta": "/.fEta",
        "phi": "/.fPhi",
        "energy": "/.fE",
        "x": "/.fX",
        "y": "/.fY",
        "z": "/.fZ",
        "j": "jerFactor",
        "o": "origIndex",
    }

    dict_ = {}
    for n in ak_array.fields:
        ak_array_ = __prepare_array(ak_array[n], n)
        if __is_rootcompat(ak_array_):
            dict_[_rename_lookup.get(n, n)] = ak_array_

    return ak.zip(dict_)


def __make_tree_maker_event_tree(events):
    """Prepare dict with simple branches and collections that uproot can write.
    
    Args:
        events (ak.Array): Events opened with the NTreeMaker schema.

    Returns:
        dict[str, ak.Array]
    """

    out = {}
    for bname in events.fields:
        if events[bname].fields:
            sub_collection = [  # Handling sub collection first
                x.replace("Counts", "")
                for x in events[bname].fields
                if x.endswith("Counts")
            ]
            if sub_collection:
                for subname in sub_collection:
                    if events[bname][subname].fields:
                        out[f"{bname}_{subname}"] = __zip_composite(
                            ak.flatten(events[bname][subname], axis=-1)
                        )
                    else:
                        out[f"{bname}_{subname}"] = ak.flatten(
                            events[bname][subname], axis=-1
                        )
            out[bname] = __zip_composite(events[bname])
        else:
            out[bname] = __prepare_array(events[bname], bname)
    return out


def write_tree_maker_root_file(output_file_name, events=None, trees={}, mode="recreate"):
    """Write events opened with the TreeMaker or NTreeMaker schema to a ROOT file.
    
    Args:
        output_file_name (str)
        events (ak.Array)
        trees (dict[str, any]): Other trees with no collections
            Keys are tree name
            Values are tree content
        mode (str): "recreate" or "update"
    """

    log.blank_line()
    log.info("Writing down output ROOT file %s" % output_file_name)
    with getattr(uproot, mode)(output_file_name) as output_file:
        if events is not None:
            output_file["Events"] = __make_tree_maker_event_tree(events)
            log.info("TTree Events saved to output file")
        for tree_name, tree in trees.items():
            output_file[tree_name] = ak.Array(tree)
            log.info("TTree %s saved to output file" % tree_name)

