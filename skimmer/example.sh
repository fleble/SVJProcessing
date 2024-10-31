
# Local copy of TreeMaker files to do fast test at PSI T3
#input_files=/work/fleble/tmp/TreeMaker_files/0_RA2AnalysisTree.root,/work/fleble/tmp/TreeMaker_files/1_RA2AnalysisTree.root

# Remote files at LPC
prefix=root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/2018/Full_11142022/PrivateSamples/SVJ_UL2018_t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8_n-1000/
input_files=${prefix}/0_RA2AnalysisTree.root,${prefix}/1_RA2AnalysisTree.root
module=analysis_configs.t_channel_pre_selection

prefix=root://t3se01.psi.ch:1094//store/t3groups/ethz-susy/darkshowers/samples/signal/PFNanoUL/UL2018/SVJleptons_mDark_scan/PFNanoAOD/
input_files_nano=${prefix}/PFNanoAOD_SVJL_mMed-1500GeV_mDark-16GeV_rinv-0.3_alpha-peak_13TeV-pythia8_part-1.root,${prefix}/PFNanoAOD_SVJL_mMed-1500GeV_mDark-16GeV_rinv-0.3_alpha-peak_13TeV-pythia8_part-2.root
module_nano=analysis_configs.s_channel_leptons_mc_pre_selection

output_file=test.root


# Example to run locally with iterative processor (one thread) on TreeMaker NTuples
python skim.py -i ${input_files} -o ${output_file} -p ${module} -y 2018 -e iterative

# Example to run locally with futures processor (multi thread) on TreeMaker NTuples
python skim.py -i ${input_files} -o ${output_file} -p ${module} -y 2018 -e futures

# Example to run with LPC HTCondor at LPC on TreeMaker NTuples
python skim.py -i ${input_files} -o ${output_file} -p ${module} -y 2018 -e dask/lpccondor -n 2 -mem 4GB

# Example to run with SLURM, e.g. at PSI T3 on TreeMaker NTuples
python skim.py -i ${input_files} -o ${output_file} -p ${module} -y 2018 -e dask/slurm -n 2 -mem 4GB -t 20:00 -q short

# Example to run locally on (PF)NanoAOD NTuples
python skim.py -i ${input_files_nano} -o ${output_file} -p ${module_nano} -y 2018 -e futures -nano  # -xsec <x> to add the cross section, replace <x> by the cross section value

# Example to run with ETP HTCondor at KIT
python skim.py -i ${input_files} -o ${output_file} -p ${module} -y 2018 -e dask/ETPCondor -n 2

