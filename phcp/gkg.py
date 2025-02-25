import logging
import shlex
import subprocess

from phcp.config import GKG_CONTAINER_PATHS


logger = logging.getLogger(__name__)


def run_gkg_command(
    cmd: list[str] | str,
    *,
    gkg_container_version,
    input_dirs: list[str] = [],
    output_dirs: list[str] = [],
    **kwargs,
) -> subprocess.CompletedProcess:
    """Run a command within a Gkg container, similarly to subprocess.run.

    Contrary to subprocess.run, the default is check=True (raise a
    CalledProcessError in case of a non-zero return code).
    """
    if "check" not in kwargs:
        kwargs["check"] = True
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    container_path = GKG_CONTAINER_PATHS[gkg_container_version]
    full_cmd = (
        [
            "singularity",
            "exec",
        ]
        + sum((["--bind", f"{dir}:{dir}:ro"] for dir in input_dirs), [])
        + sum((["--bind", f"{dir}:{dir}:rw"] for dir in output_dirs), [])
        + [container_path]
        + cmd
    )
    logger.info("Running Gkg command: %s", shlex.join(full_cmd))
    return subprocess.run(full_cmd, shell=False, **kwargs)


def run_gkg_GetMask(
    args: list[str],
    *,
    input_dirs: list[str] = [],
    output_dirs: list[str] = [],
    **kwargs,
) -> subprocess.CompletedProcess:
    return run_gkg_command(
        ["GkgExecuteCommand", "GetMask"] + args,
        gkg_container_version="2022-12-20",
        input_dirs=input_dirs,
        output_dirs=output_dirs,
        **kwargs,
    )
