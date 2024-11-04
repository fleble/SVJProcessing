#!/bin/bash

dataset_directory=/ceph/mgais/PFNanoAOD_SVJ_std2_UL2018_scouting_truth_study_inv_parts
dataset_config=dataset_configs.s_channel_scouting_signal_paths

module=analysis_configs.s_channel_scouting_pre_selection
selection_name=s_channel_scouting_pre_selection

#module=analysis_configs.t_channel_wnae_qcd_training_region
#selection_name=t_channel_wnae_qcd_training_region

#module=analysis_configs.t_channel_wnae_top_training_region
#selection_name=t_channel_wnae_top_training_region

#module=analysis_configs.t_channel_lost_lepton_control_region
#selection_name=t_channel_lost_lepton_control_region

year=2018

dataset_names=(
    #
    # Signals
    #
    mMed-700GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV

    mMed-800GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV
    mMed-800GeV_mDark-20GeV_rinv-0.5_alpha-peak_13TeV
    mMed-800GeV_mDark-20GeV_rinv-0.7_alpha-peak_13TeV

    mMed-900GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV
    mMed-900GeV_mDark-20GeV_rinv-0.5_alpha-peak_13TeV
    mMed-900GeV_mDark-20GeV_rinv-0.7_alpha-peak_13TeV

    mMed-1000GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV
    mMed-1000GeV_mDark-20GeV_rinv-0.5_alpha-peak_13TeV
    mMed-1000GeV_mDark-20GeV_rinv-0.7_alpha-peak_13TeV

    mMed-1500GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV
    mMed-1500GeV_mDark-20GeV_rinv-0.5_alpha-peak_13TeV
    mMed-1500GeV_mDark-20GeV_rinv-0.7_alpha-peak_13TeV
)

prepare_input_files_list() {

    local dataset_config=$1
    local dataset_directory=$2
    local module=$3
    local selection_name=$4
    local year=$5
    local dataset_name=$6

    echo ""
    echo "Preparing input files for dataset ${dataset_name} year ${year} and selection ${selection_name}"

    python list_dataset_files.py -d ${dataset_name} -y ${year} -c ${dataset_config} -o ${dataset_directory} -nano_scout
    python compute_unweighted_selection_efficiency.py -d ${dataset_name} -y ${year} -p ${module} -s ${selection_name} -i ${dataset_directory} -o ${dataset_directory} -n 6 -e futures -c 10000 -nano_scout
    python prepare_input_files_list.py -d ${dataset_name} -y ${year} -s ${selection_name} -i ${dataset_directory} -o ${dataset_directory} -m 50000
}


for dataset_name in ${dataset_names[@]}; do

    prepare_input_files_list ${dataset_config} ${dataset_directory} ${module} ${selection_name} ${year} ${dataset_name}

done

