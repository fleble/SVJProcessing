
# Local copy of TreeMaker files to do fast test
input_files=/work/fleble/tmp/TreeMaker_files/0_RA2AnalysisTree.root,/work/fleble/tmp/TreeMaker_files/1_RA2AnalysisTree.root

# Remote files
#input_files=root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/2018/Full_11142022/PrivateSamples/SVJ_UL2018_t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8_n-1000/0_RA2AnalysisTree.root,root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/2018/Full_11142022/PrivateSamples/SVJ_UL2018_t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8_n-1000/1_RA2AnalysisTree.root

output_file=test.root
module=analysisConfigs.example_preselection

python skim.py -i ${input_files} -o ${output_file} -m ${module} -y 2018

