#!/bin/bash

MEMORY=4GB
CORES=2
CHUNK_SIZE=10000
N_WORKERS=100
EXECUTOR=dask/lpccondor   # HTCondor at LPC
#N_WORKERS=6
#EXECUTOR=futures     # local job
FORCE_RECREATE=0   # 1 to recreate output file if it exists, 0 else
FIRST_FILE=0
LAST_FILE=-1  # Use -1 to skim all input files, or a positive integer

dataset_directory=${HOME}/nobackup/SVJ/store/datasets
output_directory=root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/signal_production_2Dscans/skims

module=analysis_configs.t_channel_pre_selection
selection_name=t_channel_pre_selection

#module=analysis_configs.t_channel_wnae_qcd_training_region
#selection_name=t_channel_wnae_qcd_training_region

#module=analysis_configs.t_channel_wnae_top_training_region
#selection_name=t_channel_wnae_top_training_region

#module=analysis_configs.t_channel_lost_lepton_control_region
#selection_name=t_channel_lost_lepton_control_region

years=(
    2016
    2016APV
    2017
    2018
)

variations=(
    nominal
#    # JEC/JER variations only for signal!
#    jec_up
#    jec_down
#    jer_up
#    jer_down
#    ue_up
#    ue_down
)

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
    TTJets_SingleLeptFromT
    TTJets_SingleLeptFromTbar
    TTJets_SingleLeptFromT_genMET-150
    TTJets_SingleLeptFromTbar_genMET-150
    TTJets_DiLept
    TTJets_DiLept_genMET-150
    #
    # WJets
    #
    WJetsToLNu_HT-400To600
    WJetsToLNu_HT-600To800
    WJetsToLNu_HT-800To1200
    WJetsToLNu_HT-1200To2500
    WJetsToLNu_HT-2500ToInf
    #
    # ZJets
    #
    ZJetsToNuNu_HT-400To600
    ZJetsToNuNu_HT-600To800
    ZJetsToNuNu_HT-800To1200
    ZJetsToNuNu_HT-1200To2500
    ZJetsToNuNu_HT-2500ToInf
    #
    # Single top
    #
    ST_s-channel_4f_hadronicDecays
    ST_s-channel_4f_leptonDecays
    ST_t-channel_antitop_5f_InclusiveDecays
    ST_t-channel_top_5f_InclusiveDecays
    ST_tW_top_5f_inclusiveDecays
    ST_tW_antitop_5f_inclusiveDecays
    #
    # Diboson
    #
    # ZZTo2Q2Nu
    # WZTo2Q2Nu
    # WZTo1L1Nu2Q
    # WWTo1L1Nu2Q
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


make_skims() {

    local dataset_directory=$1
    local module=$2
    local selection_name=$3
    local year=$4
    local variation=$5
    local dataset_name=$6
    local output_directory=$7

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
                local output_file_tmp=/uscmst1b_scratch/lpc1/3DayLifetime/${USER}/${output_file_name_tmp}

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
                    if [[ ${dataset_name} == t-channel* ]]; then
                        weight_variation_flag="--weight_variation scale pdf"
                    else
                        weight_variation_flag=""
                    fi
                    python skim.py -i ${input_files} -o ${output_file_tmp} -p ${module} -pd ${dataset_name} -y ${year} -e ${EXECUTOR} -n ${N_WORKERS} -c ${CHUNK_SIZE} --memory ${MEMORY} --cores ${CORES} -pn_tagger ${variation_flag} ${weight_variation_flag}
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
  for year in ${years[@]}; do
    for variation in ${variations[@]}; do
      make_skims ${dataset_directory} ${module} ${selection_name} ${year} ${variation} ${dataset_name} ${output_directory}
    done
  done
done

