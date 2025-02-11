import logging
import os
import shlex
import subprocess

from phcp.config import FSLDIR

logger = logging.getLogger(__name__)


def run_fsl_command(cmd: list[str] | str, **kwargs) -> subprocess.CompletedProcess:
    """Run a command within a FSL environment, similarly to subprocess.run.

    Contrary to subprocess.run, the default is check=True (raise a
    CalledProcessError in case of a non-zero return code).

    The path to FSL can be configured in the phcp.config.FSLDIR variable.
    """
    if "check" not in kwargs:
        kwargs["check"] = True
    if not isinstance(cmd, str):
        cmd = shlex.join(cmd)
    if "shell" in kwargs:
        del kwargs["shell"]
    if "env" in kwargs:
        env = kwargs["env"].copy()
        del kwargs["env"]
    else:
        env = os.environ.copy()
    env["FSLDIR"] = FSLDIR
    fsl_sh = os.path.join(FSLDIR, "etc", "fslconf", "fsl.sh")
    if not os.path.isfile(fsl_sh):
        raise RuntimeError(f"Cannot find {fsl_sh}")
    logger.info("Running FSL command: %s", cmd)
    return subprocess.run(
        f". {shlex.quote(fsl_sh)}; {cmd}", env=env, shell=True, **kwargs
    )
