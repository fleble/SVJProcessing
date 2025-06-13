from coffea.nanoevents.schemas.base import BaseSchema, zip_forms, nest_jagged_forms


class NTreeMakerSchema(BaseSchema):
    def __init__(self, base_form):
        super().__init__(base_form)
        self._remove_faulty_branches()
        self._form["contents"] = self._build_collections(self._form["contents"])

    def _remove_faulty_branches(self):
        keys = list(self._form["contents"].keys())

        for key in keys:
            # Removing subjet branches as they were affected by a bug
            if "_subjets" in key:
                # TODO: This could have been removed after April 2024 TreeMaker update...
                # To remove if making skims from scratch again!
                self._form["contents"].pop(key)

            # Removing automatically added counter branches
            branch_name = "n" + key.split("_")[0]
            if branch_name in self._form["contents"].keys():
                self._form["contents"].pop(branch_name)

        # Remove specific branches causing issues
        extra_variables_to_remove = [
            "nJetsAK8_constituentsIndex",
            "Photons_electronFakes",
            "Photons_nonPrompt",
        ]
        for var in extra_variables_to_remove:
            if var in self._form["contents"].keys():
                self._form["contents"].pop(var)

    def _build_collections(self, branch_forms):
        # Turn any special classes into the appropriate awkward form
        composite_objects = list(
            set(k.split("/")[0].rstrip("_") for k in branch_forms if "/" in k)
        )

        composite_behavior = {  # Dictionary for overriding the default behavior
            "Tracks": "LorentzVector"
        }
        for objname in composite_objects:
            components = {  # Extracting the various composit object names
                k.split(".")[-1]: k
                for k in branch_forms
                if k.startswith(objname + "/") or
                # Second case for skimming
                k.startswith(objname + "_/")
            }

            if set(components.keys()) == {
                "fPt",
                "fEta",
                "fPhi",
                "fE",
            }:
                form = zip_forms(
                    {
                        "pt": branch_forms.pop(components["fPt"]),
                        "eta": branch_forms.pop(components["fEta"]),
                        "phi": branch_forms.pop(components["fPhi"]),
                        "energy": branch_forms.pop(components["fE"]),
                    },
                    objname,
                    composite_behavior.get(objname, "PtEtaPhiELorentzVector"),
                )
                branch_forms[objname] = form
            elif set(x.split(".")[-1] for x in components) == {
                "fX",
                "fY",
                "fZ",
            }:
                form = zip_forms(
                    {
                        "x": branch_forms.pop(components["fX"]),
                        "y": branch_forms.pop(components["fY"]),
                        "z": branch_forms.pop(components["fZ"]),
                    },
                    objname,
                    composite_behavior.get(objname, "ThreeVector"),
                )
                branch_forms[objname] = form
            else:
                raise ValueError(
                    f"Unrecognized class with split branches of object {objname}: {components.values()}"
                )

        # Generating collection from branch name
        # Cannot start with "n" as this is generated from uproot counting
        collections = [k for k in branch_forms if "_" in k and not k.startswith("n")]
        collections = sorted(
            set(
                [
                    "_".join(k.split("_")[:-1])
                    for k in collections
                    if k.split("_")[-1] != "AK8"
                    # Excluding per-event variables with AK8 variants like Mjj and MT
                ]
            ),
            # Always process nested items first
            key=lambda colname: colname.count("_"),
            reverse=True,
        )

        subcollections = []

        for cname in collections:
            items = sorted(k for k in branch_forms if k.startswith(cname + "_"))
            if len(items) == 0:
                continue

            # Special pattern parsing for <collection>_<subcollection>Counts branches
            countitems = [x for x in items if x.endswith("Counts")]
            subcols = set(x[:-6] for x in countitems)  # List of subcollection names
            for subcol in subcols:
                items = [
                    k for k in items if not k.startswith(subcol) or k.endswith("Counts")
                ]
                subname = subcol[len(cname) + 1 :]
                subcollections.append(
                    {
                        "colname": cname,
                        "subcol": subcol,
                        "countname": subname + "Counts",
                        "subname": subname,
                    }
                )

            if cname not in branch_forms:
                collection = zip_forms(
                    {k[len(cname) + 1]: branch_forms.pop(k) for k in items}, cname
                )
                branch_forms[cname] = collection
            else:
                collection = branch_forms[cname]
                if not collection["class"].startswith("ListOffsetArray"):
                    raise NotImplementedError(
                        f"{cname} isn't a jagged array, not sure what to do"
                    )
                for item in items:
                    itemname = item[len(cname) + 1 :]
                    collection["content"]["contents"][itemname] = branch_forms.pop(
                        item
                    )["content"]

        for sub in subcollections:
            nest_jagged_forms(
                branch_forms[sub["colname"]],
                branch_forms.pop(sub["subcol"]),
                sub["countname"],
                sub["subname"],
            )

        return branch_forms

    @property
    def behavior(self):
        """Behaviors necessary to implement this schema"""
        from coffea.nanoevents.methods import base, vector

        behavior = {}
        behavior.update(base.behavior)
        behavior.update(vector.behavior)
        return behavior

