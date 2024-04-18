import numpy as np
import awkward as ak
from utils.particle_net import jetutils as ju
from utils.particle_net.svjgnntagger import SVJGNNTagger
from scipy.special import softmax

def run_jet_tagger(events,fjets):
    gnn_triton = SVJGNNTagger(score_tag='score',
            triton_path='triton+grpc://triton.fnal.gov:443/svj_tch_gnn/1',
            model_structure='utils.data.GNNTagger.SVJTagger',
            model_inputs='../utils/data/GNNTagger/svj.yaml',
            dec_thresh=0.999)
    # initialize model if not already done
    gnn_triton.use_triton = True
    if gnn_triton.model == None:
        gnn_triton.initialize_model()
    jets_in = ju.run_jet_constituent_matching(events,fjets)
    jets_in = ak.flatten(jets_in)
    batch_size = 1024
    nnOutput = np.array([])
    for ii in range(0,len(jets_in),batch_size):
        try:
            jets_eval = jets_in[ii:ii+batch_size]
        except:
            jets_eval = jets_in[ii:-1]
        # put data in proper format
        feature_map = gnn_triton.get_feature_map(jets_eval)
        X = gnn_triton.structure_X(jets_eval,feature_map)
        #inference with triton 
        outputs = gnn_triton.model(X)
        nnOutput = np.append(nnOutput,softmax(outputs, axis=-1)[:, 2]) # 2 is the label for SVJ_Dark, 0 = QCD, 1 = TTJets
    counts = ak.num(fjets.pt)
    return ak.unflatten(nnOutput, counts)
