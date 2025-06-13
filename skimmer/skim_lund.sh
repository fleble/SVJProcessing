inputFile=~/svj/CMSSW_10_6_29_patch1/src/TreeMaker/run/SVJ_UL2018_t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8_n-1000_RA2AnalysisTree.root
outputFile=SVJ_UL2018_t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8_n-1000_RA2AnalysisTree_SKIM.root
# inputFile=~/0_RA2AnalysisTree.root
# outputFile=svj-test.root
module=analysis_configs.t_channel_pre_selection
selection_name=t_channel_pre_selection
dataset_name=t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1
year=2018
variations=nominal
variation_flag="--variation nominal"
weight_variation_flag="--weight_variation lund"

python skim.py -i ${inputFile} -o ${outputFile} -p ${module} -pd ${dataset_name} -y ${year} -lund #-pn_tagger #${weight_variation_flag} #${variation_flag} #-e ${EXECUTOR} -n ${N_WORKERS} -c ${CHUNK_SIZE} --memory ${MEMORY} --cores ${CORES} 
