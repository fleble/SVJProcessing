import os

from coffea import processor, nanoevents
from dask_jobqueue import SLURMCluster
from distributed import Client


def __get_client(executor_name, n_workers, port=8787):
    """Return batch system client for job submission.

    Args:
        executor_name (str)
        n_workers (int)
        port (int)

    Return:
        SLURMCluster or LPCCondorCluster
    """

    job_script_prologue = [
        f'export X509_USER_PROXY={os.environ["X509_USER_PROXY"]}',
    ]

    if "slurm" in executor_name:
        cluster = SLURMCluster(
            cores=n_workers,
            processes=n_workers,
            memory="4GB",
            job_script_prologue=job_script_prologue,
        )

    elif "lpccondor" in executor_name:
        from lpcjobqueue import LPCCondorCluster

        repo_directory = os.environ["SVJ_PROCESSING_ROOT"]

        cluster = LPCCondorCluster(
            scheduler_options={
                "dashboard_address": f"{port}",
            },
            cores=1,
            memory="2GB",
            disk="100MB",
            transfer_input_files=[
                f"{repo_directory}/utils",
                f"{repo_directory}/skimmer",
                f"{repo_directory}/analysis_configs",
            ],
            log_directory=None,
            death_timeout=180,
            job_script_prologue=job_script_prologue,
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
        executor_args["client"] = __get_client(executor_name, n_workers, port)

    return executor_args

