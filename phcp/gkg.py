import logging
import os
import pathlib
import shlex
import subprocess

import numpy as np

from phcp.config import GKG_CONTAINER_PATHS


logger = logging.getLogger(__name__)


FileSpec = str | os.PathLike


def _process_container_mounts(
    input_dirs: list[FileSpec], output_dirs: list[FileSpec]
) -> tuple[list[FileSpec], list[FileSpec]]:
    # Try to cover common cases with symlinks to directories, by including both
    # the symlink and the resolved path
    input_dirs = {pathlib.Path(d).absolute() for d in input_dirs} | {
        pathlib.Path(d).resolve(strict=True) for d in input_dirs
    }
    output_dirs = {pathlib.Path(d).absolute() for d in output_dirs} | {
        pathlib.Path(d).resolve(strict=False) for d in output_dirs
    }

    # Avoid mounting the same directory twice, giving priority to parent
    # directories over their children, and to read-write mounts over read-only
    # mounts. Directories are sorted by increasing depth so we only have to go
    # through the list once, removing directories from the tail. Removing
    # elements while iterating is safe here, because they are always removed from
    # right to left.
    output_dirs = sorted(output_dirs, key=lambda p: len(p.parts))
    for i, higher_dirpath in enumerate(output_dirs):
        to_remove = []
        for j, lower_dirpath in enumerate(output_dirs[i + 1 :], start=i + 1):
            if higher_dirpath in lower_dirpath.parents:
                to_remove.append(j)
        for j in reversed(to_remove):
            del output_dirs[j]

    # Same thing with input dirs, except we begin by checking for duplicates in
    # the output dirs.
    input_dirs = sorted(input_dirs, key=lambda p: len(p.parts))
    to_remove = []
    for i, lower_dirpath in enumerate(input_dirs):
        for higher_dirpath in output_dirs:
            if (
                higher_dirpath == lower_dirpath
                or higher_dirpath in lower_dirpath.parents
            ):
                to_remove.append(i)
    for j in reversed(to_remove):
        del input_dirs[j]
    for i, higher_dirpath in enumerate(input_dirs):
        to_remove = []
        for j, lower_dirpath in enumerate(input_dirs[i + 1 :], start=i + 1):
            if higher_dirpath in lower_dirpath.parents:
                to_remove.append(j)
        for j in reversed(to_remove):
            del input_dirs[j]

    return (input_dirs, output_dirs)


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
    input_dirs, output_dirs = _process_container_mounts(input_dirs, output_dirs)
    full_cmd = (
        [
            "singularity",
            "exec",
            "--cleanenv",
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


def run_gkg_SubVolume(
    args: list[str],
    *,
    input_dirs: list[str] = [],
    output_dirs: list[str] = [],
    **kwargs,
) -> subprocess.CompletedProcess:
    return run_gkg_command(
        ["GkgExecuteCommand", "SubVolume"] + args,
        gkg_container_version="2022-12-20",
        input_dirs=input_dirs,
        output_dirs=output_dirs,
        **kwargs,
    )


def gkg_convert_gis_to_nifti(
    input_gis: str, output_nifti: str, *, verbose=False, **kwargs
) -> subprocess.CompletedProcess:
    command = [
        "GkgExecuteCommand",
        "Gis2NiftiConverter",
        "-i",
        input_gis,
        "-o",
        output_nifti,
    ]
    if verbose:
        command += ["-verbose"]

    input_dir = os.path.dirname(input_gis)
    output_dir = os.path.dirname(output_nifti)

    return run_gkg_command(
        command,
        gkg_container_version="2022-12-20",
        input_dirs=[input_dir],
        output_dirs=[output_dir],
        **kwargs,
    )


def load_GIS_image(GisFilename):
    DimFilename = GisFilename[:-3] + "dim"
    TxtInDimFile = open(DimFilename, "r").read()
    TxtFirstLine = TxtInDimFile.split("\n")[0]
    shape = tuple(np.int16(TxtFirstLine.split(" ")))  # [:-1]
    arr = np.fromfile(GisFilename, np.float32())
    arr = np.reshape(arr, shape, order="F")
    return arr


def arrange_ArrayToFlipHeaderFormat(arr):
    newarr = np.swapaxes(arr, 1, 2)
    newarr = np.flip(newarr, 1)
    newarr = np.flip(newarr, 0)
    return newarr


def dearrange_ArrayToFlipHeaderFormat(arr):
    newarr = np.flip(arr, 0)
    newarr = np.flip(newarr, 1)
    newarr = np.swapaxes(newarr, 1, 2)
    return newarr
