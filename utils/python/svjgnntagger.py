import importlib, yaml, torch

import numpy as np
import awkward as ak

from utils.data.GNNTagger.config import DataConfig
from utils.python.tritonutils import wrapped_triton
from utils.python.jetutils import *
from scipy.special import softmax


class SVJGNNTagger(object):
  """
  Class used for handling SVJ GNN Tagging using the awkward array interface
  """
  def __init__(self,
               score_tag='',
               tag='',
               title='',
               local_path=None,
               triton_path=None,
               model_structure='',
               model_inputs='',
               model=None,
               radius=0.8,
               dec_thresh=0.5,
               ipz_cut=100):
    self.score_tag = score_tag
    self.tag = tag
    self.title = title
    self.local_path = local_path
    self.triton_path = triton_path
    self.model_structure = model_structure
    self.model_inputs = model_inputs
    self.data_config = self.load_data_config()

    self.model = None
    self.device = None
    self.use_triton = None

    self.radius = radius
    self.ipz_cut = ipz_cut
    self.dec_thresh = dec_thresh

    # initialize model before inference if using local pt model and during
    # inference if using triton model
    if self.local_path:
      self.initialize_model()

  def initialize_model(self):
    """
    Initialize model for inference with local torch model or triton client
    """

    if self.triton_path:
      # check if triton server will connect
      try:
        self.model = wrapped_triton(self.triton_path)
        self.use_triton = True
      except:
        raise Exception('SVJGNNTagger error: Triton path will not load.')
    elif self.local_path:
      # check if local model will connect
      try:
        self.device = torch.device(
            "cuda" if torch.cuda.is_available() else "cpu")
        self.model = self.load_model()
        self.use_triton = False
      except:
        raise Exception('SVJGNNTagger error: Local path will not load.')
    else:
      raise Exception(
          'SVJGNNTagger error: Neither triton or local path is defined.',
          'Must define one or the other.')

    return

  def load_data_config(self):
    """
    Load in data configuration for GNN model
    """

    # get proper configuration of data
    with open(self.model_inputs, 'r') as file:
      yaml_file = yaml.safe_load(file)
    data_config = DataConfig(yaml_file, print_info=False)

    return data_config

  def load_model(self):
    """
    Load in GNN model and get it ready to evaluate
    """

    # import proper model
    mod = importlib.import_module(self.model_structure)

    # load configuration into model structure
    model, model_info = mod.get_model(self.data_config)
    model.load_state_dict(
        torch.load(self.local_path, map_location=torch.device(self.device)))
    model.eval()

    return model

  def get_feature_map(self, jets):
    """
    Put all of the possible features to be used by model in proper format
    """
    feature_map = {
        'del_eta':
        delta_eta(jets.Constituents.eta,jets.eta),
        'del_phi':
        delta_phi(jets.Constituents.phi,jets.phi),
        'log_pt':
        np.log(jets.Constituents.pt),
        'log_e':
        np.log(jets.Constituents.energy),
        'log_pt_jetpt':
        np.log(jets.Constituents.pt / jets.pt),
        'log_e_jete':
        np.log(jets.Constituents.energy / jets.energy),
        'del_r':
        delta_R(jets.Constituents.eta,jets.eta,jets.Constituents.phi,jets.phi),
    }
    # pdgId uses one-hot encoding
    pdgIds = [-13,-211,1,11,13,130,2,211,22]
    for pid in pdgIds:
        feature_map["pdgId{}".format(pid)] = np.where(jets.Constituents.PdgId == pid, 1, 0)
    feature_map['mask'] = ak.ones_like(feature_map['del_eta'])

    return feature_map

  def structure_X(self, jets, feature_map):
    """
    Put the jet data in proper format for model evaluation, X will be the
    input of the model holding all features in the proper format described
    by the data configuration file
    """
    X = {}
    for jj in range(len(self.data_config.input_names)):
      in_name = self.data_config.input_names[jj]
      in_shapes = self.data_config.input_shapes[in_name]
      X[in_name] = np.empty((len(jets), in_shapes[1], in_shapes[2]),
                            dtype=np.float32)
      for ii in range(in_shapes[1]):
        feat = self.data_config.input_dicts[in_name][ii]
        X[in_name][:, ii] = ak.to_numpy(ak.fill_none(
            ak.pad_none(feature_map[feat], in_shapes[2], axis=-1, clip=True),
            0))

      # slight modifications if using triton or local model
      if self.use_triton:
        X[f'{in_name}__{jj:d}'] = X.pop(in_name)
      else:
        X[in_name] = torch.from_numpy(X[in_name])
    return X

  def triton_evaluate(self, X):
    return self.model(X)

  def run_tag(self, events, jets_in, **kwargs):
    """
    Evaluate jets with model to get score (probability of SVJ) or binary SVJ tag
    """

    batch_size = kwargs.get('batch_size', 1024)

    # don't rerun inference if already done (when using new dec. thresh.)
    if self.score_tag in jets_in.fields:
      jets_out = jets_in[:]
      jets_out[self.tag] = jets_out[self.score_tag] > self.dec_thresh
      return jets_out

    # initialize model if not already done
    if self.model == None:
      self.initialize_model()

    # change dimensions into jet-level and get features needed
    counts = ak.num(jets_in.pt)
    jets = ak.flatten(jets_in)

    scores = np.array([])
    # split evaluation into batches
    for ii in range(0, len(jets), batch_size):

      try:
        jets_eval = jets[ii:ii + batch_size]
      except:
        jets_eval = jets[ii:-1]

      feature_map = self.get_feature_map(jets_eval)
      X = self.structure_X(jets_eval, feature_map)

      # inference slightly different depending on model using
      if self.use_triton:
        outputs = self.triton_evaluate(X)
        scores = np.append(scores, softmax(outputs, axis=-1)[:, 0])
      else:
        inputs = [X[k].to(self.device) for k in self.data_config.input_names]
        with torch.no_grad():
          outputs = self.model(*inputs)
        scores = np.append(scores,
                           torch.softmax(outputs, dim=1)[:, 0].detach().numpy())
    # change scores to event-based feature
    scores = ak.unflatten(scores, counts)

    jets_out = jets_in[:]
    jets_out[self.score_tag] = scores
    jets_out[self.tag] = scores > self.dec_thresh
    return jets_out

  @property
  def is_dr4(self):
    """
    Returning whether this tagger is a DR8 tagging scheme or not
    """
    return self.radius == 0.4

  @property
  def acc_cutname(self):
    """
    Creating a string to identify the jet-track association cut values.
    """
    r = 4 if self.is_dr4 else 8
    return f'DR{r}_ipz{self.ipz_cut:.1f}'.replace('.', 'p')

  @property
  def acc_cuttitle(self):
    """
    Creating a plot ready string to be displayed in the jet-track association
    cut-values.
    """
    return f'$\Delta R$<{self.radius}, $IP_{{z}}$<{self.ipz_cut} cm'
