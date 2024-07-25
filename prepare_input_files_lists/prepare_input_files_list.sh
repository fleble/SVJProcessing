#!/bin/bash

dataset_directory=/work/cazzanig/datasets
dataset_config=dataset_configs.s_channel_leptons_datasets_paths

module=analysis_configs.s_channel_leptons_mc_pre_selection
selection_name=s_channel_leptons_pre_selection

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
    
    #
    # Backgrounds
    #
    # QCD
    #
    #qcd_pt-470to600GeV
    #qcd_pt-600to800GeV
    #qcd_pt-800to1000GeV
    #qcd_pt-1000to1400GeV
    #qcd_pt-1400to1800GeV
    #qcd_pt-1800to2400GeV
    #qcd_pt-2400to3200GeV
    #qcd_pt-3200toInfGeV
    #
    # TTJets
    #
    #TTJets_ht-600to800GeV
    #TTJets_ht-800to1200GeV
    #TTJets_ht-1200to2500GeV
    #TTJets_ht-2500toInfGeV
    #
    # WJets
    #
    WJetsToLNu_ht-100to200GeV
    WJetsToLNu_ht-200to400GeV
    WJetsToLNu_ht-400to600GeV
    WJetsToLNu_ht-600to800GeV
    #WJetsToLNu_ht-800to1200GeV
    #WJetsToLNu_ht-1200to2500GeV
    #WJetsToLNu_ht-2500toInfGeV
    #
    # ZJets
    #
    
    #
    # Single top
    #
    
    #
    # Diboson
    #
    
    #
    ##################################
    # DATA
    ##################################
    #
    
    JetHT
    MET
    HTMHT
    EGamma
    SingleElectron
    SingleMuon

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

    #python list_dataset_files.py -d ${dataset_name} -y ${year} -c ${dataset_config} -o ${dataset_directory} -nano
    python compute_unweighted_selection_efficiency.py -d ${dataset_name} -y ${year} -p ${module} -s ${selection_name} -i ${dataset_directory} -o ${dataset_directory} -n 6 -e futures -c 10000 -nano
    python prepare_input_files_list.py -d ${dataset_name} -y ${year} -s ${selection_name} -i ${dataset_directory} -o ${dataset_directory} -m 10000 
}


for dataset_name in ${dataset_names[@]}; do

    prepare_input_files_list ${dataset_config} ${dataset_directory} ${module} ${selection_name} ${year} ${dataset_name}

done

