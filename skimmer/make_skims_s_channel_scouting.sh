#!/bin/bash

MEMORY=4GB
CORES=1
CHUNK_SIZE=10000
N_WORKERS=20
EXECUTOR=dask/ETPCondor   # HTCondor at KIT ETP
#N_WORKERS=6
#EXECUTOR=futures     # local job
FORCE_RECREATE=0   # 1 to recreate output file if it exists, 0 else
FIRST_FILE=0
LAST_FILE=-1  # Use -1 to skim all input files

dataset_directory=/ceph/mgais/PFNanoAOD_SVJ_std2_UL2018_scouting_truth_study_inv_parts

module=analysis_configs.s_channel_scouting_pre_selection
selection_name=s_channel_scouting_pre_selection

#module=analysis_configs.t_channel_wnae_qcd_training_region
#selection_name=t_channel_wnae_qcd_training_region

#module=analysis_configs.t_channel_wnae_top_training_region
#selection_name=t_channel_wnae_top_training_region

#module=analysis_configs.t_channel_lost_lepton_control_region
#selection_name=t_channel_lost_lepton_control_region

year=2018

variation=nominal  # nominal jec_up jec_down jer_up jer_down

# Output directory for nominal samples - no variation of the uncertainties
#output_directory=root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/${year}/Full/PrivateSkims/${variation}
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


make_skims() {

    local dataset_directory=$1
    local module=$2
    local selection_name=$3
    local year=$4
    local dataset_name=$5
    local output_directory=$6

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
                local output_file_tmp=/tmp/${USER}/${output_file_name_tmp}

                echo ""
                echo "Making skim file ${output_file}"

                local output_redirector=$(echo ${output_file} | cut -d/ -f 1-4)
                local output_file_path=$(echo ${output_file} | cut -d/ -f 4-)
                xrdfs ${output_redirector} ls ${output_file_path} > /dev/null 2>&1
                if [ "$?" != "0" ] || [ "${FORCE_RECREATE}" == "1" ]; then
	            if [ "${variation}" == "nominal" ]; then
		        variation_flag=''
		    else
		        variation_flag="--variation ${variation}"
		    fi
                    python skim.py -i ${input_files} -o ${output_file_tmp} -p ${module} -pd ${dataset_name} -y ${year} -nano_scout -e ${EXECUTOR} -n ${N_WORKERS} -c ${CHUNK_SIZE} --memory ${MEMORY} --cores ${CORES} -pn_tagger ${variation_flag[@]}
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

    make_skims ${dataset_directory} ${module} ${selection_name} ${year} ${dataset_name} ${output_directory}

done

