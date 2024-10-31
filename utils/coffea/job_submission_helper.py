import os

from coffea import processor, nanoevents
from dask_jobqueue import SLURMCluster, HTCondorCluster
from distributed import Client


def __get_client(executor_name, n_workers, cores, memory, disk, time, partition, port=8787):
    """Return batch system client for job submission.

    Args:
        executor_name (str)
        n_workers (int)
        cores (int)
        memory (str)
        disk (str)
        port (int)

    Return:
        SLURMCluster or LPCCondorCluster
    """

    job_script_prologue = [
        f'export X509_USER_PROXY={os.environ["X509_USER_PROXY"]}',
    ]

    if "slurm" in executor_name:
        cluster = SLURMCluster(
            cores=cores,
            memory=memory,
            log_directory=f"/work/{os.environ['USER']}/tmp/logs",
            job_script_prologue=job_script_prologue,
            walltime=time,
            queue=partition,
        )

    elif "lpccondor" in executor_name:
        from lpcjobqueue import LPCCondorCluster

        repo_directory = os.environ["SVJ_PROCESSING_ROOT"]

        cluster = LPCCondorCluster(
            scheduler_options={
                "dashboard_address": f"{port}",
            },
            cores=cores,
            memory=memory,
            disk=disk,
            transfer_input_files=[
                f"{repo_directory}/utils",
                f"{repo_directory}/skimmer",
                f"{repo_directory}/analysis_configs",
            ],
            log_directory=f"/uscmst1b_scratch/lpc1/3DayLifetime/{os.environ['USER']}/logs",
            death_timeout=180,
            job_script_prologue=job_script_prologue,
        )

    elif "ETP" in executor_name:
        repo_directory = os.environ["SVJ_PROCESSING_ROOT"]
        
        cluster = HTCondorCluster(
            scheduler_options={
                "port": "3719",
                #"dashboard_address": f"3719",
            },
            cores=cores,
            memory=memory,
            disk=disk,
            log_directory=f"{repo_directory}/logs",
            death_timeout=180,
            job_script_prologue=job_script_prologue,
            job_extra_directives = {
                "Universe":"docker", 
                "docker_image":"mschnepf/slc7-condocker", 
                "docker_network_type":"host", 
                "transfer_input_files": f"{repo_directory}/utils,{repo_directory}/skimmer,{repo_directory}/analysis_configs",
                "accounting_group":"cms.higgs"
            }
        )

        
    cluster.scale(jobs=n_workers)
    #cluster.adapt(minimum=6, maximum=250)
    client = Client(cluster)
    client.wait_for_workers(1)

    return client


def get_executor(executor_name):
    """Return coffea executor corresponding to the input name.

    Args:
        executor_name (str): Allowed names are:
           iterative
           futures
           dask/slurm
           dask/condor
           dask/ETPCondor

    Returns:
        coffea.processor
    """

    if executor_name in ["futures", "iterative"]:
        if executor_name == "iterative":
            executor = processor.iterative_executor
        else:
            executor = processor.futures_executor

    elif "dask" in executor_name:
        executor = processor.dask_executor

    return executor


def get_executor_args(
        executor_name,
        n_workers,
        cores,
        memory,
        disk,
        time, 
        partition,
        schema_name=None,
        skip_bad_files=True,
        port=8787,
    ):
    """Make executor args dictionary.

    Args:
        executor_name (str): Allowed names are:
           iterative
           futures
           dask/slurm
           dask/lpccondor
           dask/ETPCondor
        schema_name (str): e.g. BaseSchema
        n_workers (int)
        skip_bad_files (bool)

    Return
        dict
    """

    executor_args = {
        "skipbadfiles": skip_bad_files,
    }

    if schema_name is not None:
        executor_args["schema"] = getattr(nanoevents, schema_name)

    if executor_name in ["futures", "iterative"]:
        executor_args["workers"] = n_workers

    else:
        executor_args["retries"] = 5
        executor_args["client"] = __get_client(
            executor_name, n_workers, cores, memory, disk,time ,partition , port,
        )

    return executor_args

