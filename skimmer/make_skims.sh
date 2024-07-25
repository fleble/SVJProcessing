#!/bin/bash

MEMORY=6GB
TIME=12:00:00
PARTITION=standard
CORES=2
CHUNK_SIZE=10000
N_WORKERS=150
#EXECUTOR=dask/lpccondor    # HTCondor at LPC
EXECUTOR=dask/slurm     # local job
#EXECUTOR=futures     # local job
#N_WORKERS=6
FORCE_RECREATE=0   # 1 to recreate output file if it exists, 0 else
FIRST_FILE=0
LAST_FILE=-1  # Use -1 to skim all input files

dataset_directory=/work/cazzanig/datasets

module=analysis_configs.s_channel_leptons_mc_pre_selection
selection_name=s_channel_leptons_pre_selection

#module=analysis_configs.t_channel_wnae_qcd_training_region
#selection_name=t_channel_wnae_qcd_training_region

#module=analysis_configs.t_channel_wnae_top_training_region
#selection_name=t_channel_wnae_top_training_region
 
#module=analysis_configs.t_channel_lost_lepton_control_region
#selection_name=t_channel_lost_lepton_control_region

year=2018

# Output directory for nominal samples - no variation of the uncertainties
#output_directory=root://t3dcachedb03.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/cazzanig/schannel_leptons_skims/QCD/
output_directory=root://t3dcachedb03.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/cazzanig/schannel_leptons_skims/WJets/


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
    #WJetsToLNu_ht-100to200GeV
    #WJetsToLNu_ht-200to400GeV
    #WJetsToLNu_ht-400to600GeV
    #WJetsToLNu_ht-600to800GeV
    WJetsToLNu_ht-800to1200GeV
    WJetsToLNu_ht-1200to2500GeV
    WJetsToLNu_ht-2500toInfGeV
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
    #JetHT
)

cross_sections=(
    #
    # Signals
    #
    
    #
    # Backgrounds
    #
    # QCD
    #
    #551.60
    #156.4
    #26.24
    #7.477
    #0.642
    #0.08671
    #0.00523
    #0.000132
    #
    # TTJets
    #
    #2.423
    #0.9818
    #0.1714
    #0.001966
    #
    # WJets
    #
    #
    #
    #
    #
    4.933
    1.16 
    0.02627
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
    #JetHT
)

make_skims() {

    local dataset_directory=$1
    local module=$2
    local selection_name=$3
    local year=$4
    local dataset_name=$5
    local output_directory=$6
    local xsec=$7

    # Path automatically built when preparing input files lists
    local suffix_dir=${year}/${selection_name}/${dataset_name}
    local files_list_directory=${dataset_directory}/skim_input_files_list/${suffix_dir}
    local output_directory=${output_directory}/${suffix_dir}

    local output_redirector=$(echo ${output_directory} | cut -d/ -f 1-4)
    local output_dir=$(echo ${output_directory} | cut -d/ -f 4-)
    xrdfs ${output_redirector} ls ${output_dir} > /dev/null 2>&1
    if [ "$?" != "0" ]; then
        xrdfs ${output_redirector} mkdir -p ${output_dir}
    fi

    i_file=-1
    for files_list in $(ls ${files_list_directory} | sort -V); do
        ((i_file++))
        if [ ${i_file} -le ${LAST_FILE} ] || [ "${LAST_FILE}" == "-1" ]; then
            if [ ${i_file} -ge ${FIRST_FILE} ]; then

                local input_files=${files_list_directory}/${files_list}
                local output_file=${output_directory}/${files_list/.txt/.root}
                local output_file_name_tmp=$(echo ${ouput_file}_$(date +"%Y%m%d-%H%M%S") | shasum | cut -d " " -f1).root
                local output_file_tmp=/work/${USER}/tmp/${output_file_name_tmp}   

                echo ""
                echo "Making skim file ${output_file}"

                local output_redirector=$(echo ${output_file} | cut -d/ -f 1-4)
                local output_file_path=$(echo ${output_file} | cut -d/ -f 4-)
                xrdfs ${output_redirector} ls ${output_file_path} > /dev/null 2>&1
                if [ "$?" != "0" ] || [ "${FORCE_RECREATE}" == "1" ]; then
                    python skim.py -i ${input_files} -o ${output_file_tmp} -p ${module} -pd ${dataset_name} -y ${year} -e ${EXECUTOR} -n ${N_WORKERS} -c ${CHUNK_SIZE} --memory ${MEMORY} --cores ${CORES} --walltime ${TIME} --queue ${PARTITION} -xsec ${xsec}   -nano
                    xrdcp -f ${output_file_tmp} ${output_file}
                    echo ${output_file} has been saved.
                    rm ${output_file_tmp}
                else
                    echo ${output_file} already exists and FORCE_RECREATE is 0. Skipping.
                fi
            fi
        fi
    done
}


for dataset_name in ${dataset_names[@]}; do
    make_skims ${dataset_directory} ${module} ${selection_name} ${year} ${dataset_name} ${output_directory} ${cross_sections[${i}]}
    ((i++))

done

