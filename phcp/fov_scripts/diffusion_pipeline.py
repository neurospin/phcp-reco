import logging
import sys

import phcp.gkg

logger = logging.getLogger(__name__)


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        "Thin wrapper around the Gkg post-mortem diffusion pipeline."
    )

    parser.add_argument(
        "-s",
        "--subjectJsonFileName",
        required=True,
        help="Subject json dictionary filename, e.g. "
        "derivatives/gkg-Pipeline/GkgPipelineDescriptions/sub-${sub}.json",
    )
    parser.add_argument(
        "-t",
        "--taskJsonFileName",
        required=True,
        help="Tasks json dictionary filename. See the repository README for more details",
    )
    parser.add_argument(
        "-g",
        "--gkgpipelineJsonFilename",
        required=True,
        help="GkgPipelineDescription json dictionary filename",
    )
    parser.add_argument(
        "-o",
        "--outputDirectory",
        required=True,
        help="Output directory of the pipeline, e.g. fov/derivatives/gkg-Pipeline",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Show as much information as possible",
    )

    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line(argv)
    return phcp.gkg.gkg_run_diffusion_pipeline(
        args.subjectJsonFileName,
        args.taskJsonFileName,
        args.gkgpipelineJsonFilename,
        args.outputDirectory,
        verbose=args.verbose,
    ).returncode


if __name__ == "__main__":
    sys.exit(main())
