GKG_CMDLINE = [
    "singularity",
    "exec",
    "--bind",
    "/neurospin:/neurospin:rw",
    "/neurospin/phcp/code/gkg/2022-12-20_gkg/2022-12-20_gkg.sif",
    # "/neurospin/phcp/code/gkg/2024-10-26_gkg/2024-10-26_gkg.sif",
]

FSLDIR = "/drf/local/fsl-6.0.7.13"
