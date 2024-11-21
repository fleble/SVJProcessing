#!/bin/bash

MEMORY=10GB
TIME=12:00:00
PARTITION=standard
CORES=2
CHUNK_SIZE=1000
N_WORKERS=150
#EXECUTOR=dask/lpccondor    # HTCondor at LPC
EXECUTOR=futures         # run locally
FORCE_RECREATE=0   # 1 to recreate output file if it exists, 0 else
FIRST_FILE=0
LAST_FILE=-1  # Use -1 to skim all input files

dataset_directory=/work/cazzanig/datasets_svj_leptons_presel_testSYS/

module=analysis_configs.s_channel_leptons_mc_pre_selection
selection_name=s_channel_leptons_pre_selection_sys

#do not specify any path, just file name in directory data
pfnano_corrections_file=/work/cazzanig/SVJProcessing_sys/SVJProcessing/data/corrections_2024-11-19_11-13-38_jme_corr.coffea

years=(
    #2016
    #2016APV
    #2017
    2018
)

output_directory=root://t3dcachedb03.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/cazzanig/schannel_leptons_skims_nohem_mc_2018_sys_test/SVJleptons_signal/


dataset_names=(
    #
    # Signals
    #
    mMed-2000GeV_mDark-8GeV_rinv-0.3
    #
)

add_weights_variations=1

variations=(
    #nominal
#    # JEC/JER variations only for signal!
#    jec_up
#    jec_down
#    jer_up
#    jer_down
     unclEn_up
)


cross_sections=(
    #
    # Signals
    #
    1.875
)


make_skims() {

    local dataset_directory=$1
    local module=$2
    local selection_name=$3
    local year=$4
    local variation=$5
    local dataset_name=$6
    local output_directory=$7
    local xsec=$8

    # Path automatically built when preparing input files lists
    local files_list_directory=${dataset_directory}/skim_input_files_list/${year}/${selection_name}/${dataset_name}
    local output_directory=${output_directory}/${year}/${selection_name}/${variation}/${dataset_name}

    if [[ "${output_directory}"  == "root://"* ]]; then  # if the output has a redirector
        local output_redirector=$(echo ${output_directory} | cut -d/ -f 1-4)
        local output_dir=$(echo ${output_directory} | cut -d/ -f 4-)
        xrdfs ${output_redirector} ls ${output_dir} > /dev/null 2>&1
        if [ "$?" != "0" ]; then
            xrdfs ${output_redirector} mkdir -p ${output_dir}
        fi
    else
        local output_redirector="none"
        if [ ! -d ${output_directory} ]; then
            mkdir -p ${output_directory}
        fi
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

                if [ "${output_redirector}" == "none" ]; then
                    ls ${output_file} > /dev/null 2>&1
                else
                    local output_file_path=$(echo ${output_file} | cut -d/ -f 4-)
                    xrdfs ${output_redirector} ls ${output_file_path} > /dev/null 2>&1
                fi
                if [ "$?" != "0" ] || [ "${FORCE_RECREATE}" == "1" ]; then
                    if [ "${variation}" == "nominal" ]; then
                        variation_flag=""
                    else
                        variation_flag="--variation ${variation}"
                    fi
                    #if [[ ${dataset_name} == t-channel* ]]; then
                    if [ ${add_weights_variations} == 1 ]; then
                        weight_variation_flag="--weight_variation scale pdf"
                    else
                        weight_variation_flag=""
                    fi
                    python skim.py -i ${input_files} -o ${output_file_tmp} -p ${module} -pd ${dataset_name} -y ${year} -e ${EXECUTOR} -n ${N_WORKERS} -c ${CHUNK_SIZE} --memory ${MEMORY} --cores ${CORES} -pn_tagger ${variation_flag} ${weight_variation_flag} -xsec ${xsec} -corrfile ${pfnano_corrections_file} -nano -mc  
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


n_datasets=${#dataset_names[@]}

for ((i=0; i<$n_datasets; i++)); do
    dataset_name=${dataset_names[i]}
    cross_section=${cross_sections[i]}
  for year in ${years[@]}; do
    for variation in ${variations[@]}; do
      make_skims ${dataset_directory} ${module} ${selection_name} ${year} ${variation} ${dataset_name} ${output_directory} ${cross_section}
    done
  done
done