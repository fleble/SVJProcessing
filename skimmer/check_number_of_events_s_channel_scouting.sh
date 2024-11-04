#!/bin/bash

N_WORKERS=6

dataset_directory=/ceph/mgais/PFNanoAOD_SVJ_std2_UL2018_scouting_truth_study_inv_parts

selection_name=s_channel_scouting_pre_selection
#selection_name=t_channel_wnae_qcd_training_region
#selection_name=t_channel_wnae_top_training_region
#selection_name=t_channel_lost_lepton_control_region

year=2018

# Output directory for nominal samples - no variation of the uncertainties
output_directory=root://cmsdcache-kit-disk.gridka.de:1094//store/user/mgaisdor/Condor_skimmed_SVJ_std2_UL2018_scouting_truth_study_inv_parts/


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

check_number_of_events() {

    local dataset_directory=$1
    local selection_name=$2
    local year=$3
    local dataset_name=$4
    local output_directory=$5

    # Path automatically built when preparing input files lists
    local files_list=${dataset_directory}/files_list/${year}/${dataset_name}.csv
    local output_directory=${output_directory}/${year}/${selection_name}/${dataset_name}

    python check_number_of_events.py -i ${files_list} -o ${output_directory} -n ${N_WORKERS} -nano
}


for dataset_name in ${dataset_names[@]}; do

    check_number_of_events ${dataset_directory} ${selection_name} ${year} ${dataset_name} ${output_directory}

done

