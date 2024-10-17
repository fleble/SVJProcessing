import os

from tests import helper


def __prepare_input_files_lists(
        dataset_name,
        dataset_directory,
        dataset_config,
        year,
        analysis_config,
        selection_name,
    ):

    n_workers = 4
    chunk_size = 10000
    executor = "futures"

    bash_command = f"python {os.environ['SVJ_PROCESSING_ROOT']}/prepare_input_files_lists/list_dataset_files.py -d {dataset_name} -y {year} -c {dataset_config} -o {dataset_directory} -nano_scout"
    helper.test_command(bash_command)

    bash_command = f"python {os.environ['SVJ_PROCESSING_ROOT']}/prepare_input_files_lists/compute_unweighted_selection_efficiency.py -d {dataset_name} -y {year} -p {analysis_config} -s {selection_name} -i {dataset_directory} -o {dataset_directory} -e {executor} -n {n_workers} -c {chunk_size} -nano_scout"
    helper.test_command(bash_command)

    bash_command = f"python {os.environ['SVJ_PROCESSING_ROOT']}/prepare_input_files_lists/prepare_input_files_list.py -d {dataset_name} -y {year} -s {selection_name} -i {dataset_directory} -o {dataset_directory} -m 5000"
    helper.test_command(bash_command)


def test_execution():

    years = ["2018"]
    dataset_name = "mMed-1500GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV"
    dataset_config = "dataset_configs.tests_schannel_scouting"
    analysis_config = "analysis_configs.s_channel_scouting_pre_selection"
    selection_name = "s_channel_scouting_pre_selection"

    dataset_directory = helper.get_temporary_directory()
    sub_directory = helper.run_bash_command(f'echo $(date +"%Y%m%d-%H%M%S") | shasum | cut -d " " -f1')
    dataset_directory += f"/{sub_directory}"

    for year in years:
        __prepare_input_files_lists(
            dataset_name,
            dataset_directory,
            dataset_config,
            year,
            analysis_config,
            selection_name,
        )

        assert os.path.exists(f"{dataset_directory}/files_list/{year}/{dataset_name}.csv")
        assert os.path.exists(f"{dataset_directory}/selections/{year}/{selection_name}/{dataset_name}.txt")
        assert os.path.exists(f"{dataset_directory}/skim_input_files_list/{year}/{selection_name}/{dataset_name}/part-0.txt")


test_execution()

