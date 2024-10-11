import re


def pdg_id(particle):

    pdg_id_dict = {
	## Quarks
	"d": 1,
	"u": 2,
	"s": 3,
	"c": 4,
	"b": 5,
	"t": 6,

	## Leptons
	"e-"    : 11,
	"e+"    : -11,
	"mu-"   : 13,
	"mu+"   : -13,
	"tau-"  : 15,
	"tau+"  : -15,
	"nu_e"  : 12,
	"nu_mu" : 14,
	"nu_tau": 16,

	## Gauge and Higgs bosons
	"g"     : 21,
	"photon": 22,
	"Z"     : 23,
	"W+"    : 24,
	"W-"    : -24,
	"H"     : 25,

	## Light I=1 mesons
	"pi0" : 111,
	"pi+" : 211,
	"pi-" : -211,
	"rho0": 113,
	"rho+": 213,
	"rho-": -213,

	## Strange mesons
	"KL0": 130,
    }

    d = pdg_id_dict   # shorthand

    pdg_id_dict["quarks"] = [d["d"], d["u"], d["s"], d["c"], d["b"], d["t"]]

    pdg_id_dict["electron"] = [d["e-"], d["e+"]]
    pdg_id_dict["muon"] = [d["mu-"], d["mu+"]]
    pdg_id_dict["charged_leptons"] = [d["e-"], d["e+"], d["mu-"], d["mu+"], d["tau-"], d["tau+"]]
    pdg_id_dict["neutral_leptons"] = [d["nu_e"], d["nu_mu"], d["nu_tau"]]
    pdg_id_dict["leptons"] =  pdg_id_dict["charged_leptons"] + pdg_id_dict["neutral_leptons"]

    pdg_id_dict["light_mesons"] = [d["pi0"], d["pi+"], d["pi-"], d["rho0"], d["rho+"], d["rho-"]]
    pdg_id_dict["strange_mesons"] = [d["KL0"]]
    pdg_id_dict["mesons"] = d["light_mesons"] + d["strange_mesons"]

    pdg_id_dict["hadrons"] = d["mesons"]


    return pdg_id_dict[particle]


def jet_name_to_jet_radius(jet_name):
    """Return the jet radius corresponding to a jet name.

    Either radius found from the dict
    or searching for a number in the string and dividing it by 10.

    Args:
        jet_name (str): e.g. AK8, FatJet, AK8Jet, ...

    Returns:
        float
    """
    
    jet_name_to_jet_radius_dict = {
        "ak4": 0.4,
        "ak4jet": 0.4,
        "jet": 0.4,
        "ak8": 0.8,
        "ak8jet": 0.8,
        "fatjet": 0.8,
    }

    if jet_name.lower() in jet_name_to_jet_radius_dict.keys():
        radius = jet_name_to_jet_radius_dict[jet_name.lower()]

    else:
        if bool(re.search(r'\d', jet_name)):
            radius = float(re.sub('\D', '', jet_name))/10

    return radius