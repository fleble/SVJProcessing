#!/bin/bash

MEMORY=4GB
CORES=2
CHUNK_SIZE=10000
N_WORKERS=50
EXECUTOR=dask/lpccondor   # HTCondor at LPC
#N_WORKERS=6
#EXECUTOR=futures     # local job
FORCE_RECREATE=0   # 1 to recreate output file if it exists, 0 else

dataset_directory=${HOME}/nobackup/SVJ/store/datasets
module=analysis_configs.t_channel_pre_selection
selection_name=t_channel_pre_selection
year=2018
output_directory=root://cmseos.fnal.gov//store/user/${USER}/SVJ/tChannel/analysisSkims


dataset_names=(
    #
    # Signals
    #
    t-channel_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1
    t-channel_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1
    t-channel_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1
    t-channel_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1
    t-channel_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1
    t-channel_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1
    t-channel_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1
    t-channel_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1
    t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1
    t-channel_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1
    t-channel_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1
    t-channel_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1
    t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1
    #
    # Backgrounds
    #
    # QCD
    #
    QCD_Pt_170to300
    QCD_Pt_300to470
    QCD_Pt_470to600
    QCD_Pt_600to800
    QCD_Pt_800to1000
    QCD_Pt_1000to1400
    QCD_Pt_1400to1800
    QCD_Pt_1800to2400
    QCD_Pt_2400to3200
    QCD_Pt_3200toInf
    #
    # TTJets
    #
    TTJets
    TTJets_HT-600to800
    TTJets_HT-800to1200
    TTJets_HT-1200to2500
    TTJets_HT-2500toInf
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
        if [ ${i_file} -le ${LAST_FILE} ]; then
            if [ ${i_file} -ge ${FIRST_FILE} ]; then

                local input_files=${files_list_directory}/${files_list}
                local output_file=${output_directory}/${files_list/.txt/.root}
                local output_file_name_tmp=$(echo ${ouput_file}_$(date +"%Y%m%d-%H%M%S") | shasum | cut -d " " -f1).root
                local output_file_tmp=/uscmst1b_scratch/lpc1/3DayLifetime/${USER}/${output_file_name_tmp}

                echo ""
                echo "Making skim file ${output_file}"

                local output_redirector=$(echo ${output_file} | cut -d/ -f 1-4)
                local output_file_path=$(echo ${output_file} | cut -d/ -f 4-)
                xrdfs ${output_redirector} ls ${output_file_path} > /dev/null 2>&1
                if [ "$?" != "0" ] || [ "${FORCE_RECREATE}" == "1" ]; then
                    python skim.py -i ${input_files} -o ${output_file_tmp} -p ${module} -y ${year} -e ${EXECUTOR} -n ${N_WORKERS} -c ${CHUNK_SIZE} --memory ${MEMORY} --cores ${CORES}
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

