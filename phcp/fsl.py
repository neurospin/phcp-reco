import logging
import os
import shlex
import subprocess

from phcp.config import FSLDIR

logger = logging.getLogger(__name__)


def run_fsl_command(cmd: list[str] | str, **kwargs):
    """Run a command within a FSL environment, similarly to subprocess.run"""
    if not isinstance(cmd, str):
        cmd = shlex.join(cmd)
    env = dict(kwargs.get("env", os.environ))
    env["FSLDIR"] = FSLDIR
    fsl_sh = os.path.join(FSLDIR, "etc", "fslconf", "fsl.sh")
    if not os.path.isfile(fsl_sh):
        raise RuntimeError(f"Cannot find {fsl_sh}")
    logger.info("Running FSL command: %s", cmd)
    return subprocess.run(
        f". {shlex.quote(fsl_sh)}; {cmd}", env=env, shell=True, **kwargs
    )
