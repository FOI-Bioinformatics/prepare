# PREPARE documentation

Index of the documentation in this folder. The project overview, installation and the command table are in the [main README](../README.md).

## User guides

- [pipeline_overview.md](pipeline_overview.md): the pipeline stages from screening to assembly, the commands of each stage, Branchwater input and the environment commands.
- [reads.md](reads.md): choosing datasets, `sra_info`, `download_sra` with its completeness verdicts, `sra_validate`, `sra_profile` and `sra_report`.
- [project_state.md](project_state.md): the registry, `status`, `results_table`, `blacklist`, `runs` and `project_report`.
- [shared_store.md](shared_store.md): the shared data store and the `store_*` commands.
- [extraction_assembly.md](extraction_assembly.md): `extract_target_reads`, assembly options and assembly markers.
- [analysis_commands.md](analysis_commands.md): diversity, taxonomy, exploration and plotting commands with their input formats.
- [configuration.md](configuration.md): every runtime setting with its default, environment variable and flag, plus logging and the `doctor` checks.
- [hpc.md](hpc.md): SLURM array jobs, walltime signals, exit codes, resubmission, and a shared store on NFS.

## Developer guides

- [ARCHITECTURE.md](ARCHITECTURE.md): layers, design decisions, locking, exit codes and logging policy.
- [packaging.md](packaging.md): building, releasing and distributing the package.
- [../CONTRIBUTING.md](../CONTRIBUTING.md): development workflow and quality checks.
- [../CHANGELOG.md](../CHANGELOG.md): release history.
