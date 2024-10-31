
# Local copy of TreeMaker files to do fast test at PSI T3
#input_files=/work/fleble/tmp/TreeMaker_files/0_RA2AnalysisTree.root,/work/fleble/tmp/TreeMaker_files/1_RA2AnalysisTree.root

# Remote files at LPC
prefix=root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/2018/Full_11142022/PrivateSamples/SVJ_UL2018_t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8_n-1000/
input_files=${prefix}/0_RA2AnalysisTree.root,${prefix}/1_RA2AnalysisTree.root

output_file=test.root
module=analysis_configs.t_channel_pre_selection


# Example to run locally with iterative processor
python skim.py -i ${input_files} -o ${output_file} -p ${module} -y 2018 -e iterative

# Example to run with LPC HTCondor at LPC
python skim.py -i ${input_files} -o ${output_file} -p ${module} -y 2018 -e dask/lpccondor -n 2

# Example to run with SLURM, e.g. at PSI T3
python skim.py -i ${input_files} -o ${output_file} -p ${module} -y 2018 -e dask/slurm -n 2

# Example to run with ETP HTCondor at KIT
python skim.py -i ${input_files} -o ${output_file} -p ${module} -y 2018 -e dask/ETPCondor -n 2

