# Template Config Files

## jobflow_config
If you wish to use Jobflow, refer to the [Jobflow docs](https://materialsproject.github.io/jobflow/jobflow.settings.html?highlight=jobflow_config_file#jobflow.settings.JobflowSettings) for how to make a `jobflow.yaml` file to tell Jobflow where to store your calculation results. A sample `jobflow.yaml` file is included in the `jobflow_config` folder. Replace all the `???` entries with your own database information.

## fw_config
If you wish to use Fireworks, you'll need to fill out the files in `fw_config`, updating the various paths and replacing the `???` entries with your own information. This mainly includes:
- In `FW_config.yaml`, update `CONFIG_FILE_DIR` to point to the directory where you have your `FW_config.yaml` file.
- In `my_launchpad.yaml`, you will need to store your database credentials like you did for `jobflow.yaml`.
- In `my_qadapter.yaml`, you will need to provide the relevant information on how to submit a job to your cluster. An example is provided for SLURM. For VASP, note that you need to export the `VASP_PARALLEL_CMD` variable in your pre-rocket, as noted above.

For additional guidance, refer to the FireWorks setup instructions with [Atomate2](https://materialsproject.github.io/atomate2/user/fireworks.html#atomate2-fireworks).

## quacc_config
### VASP
The `vasp_custodian_settings.yaml` file found [here](https://github.com/arosen93/quacc/blob/main/quacc/defaults/custodian_settings/vasp_custodian_settings.yaml) contains the default settings for how QuAcc should run Custodian-powered VASP. If you wish to change them, you can copy this `.yaml` file and specify its path in the `VASP_CUSTODIAN_SETTINGS` environment variable. Generally though, you do not need to make any changes to this file other than `vasp_cmd` and `vasp_gamma_cmd` if your executables are named differently than `vasp_std` and `vasp_gam`.

Note that `VASP_PARALLEL_CMD` in the `vasp_custodian_settings.yaml` file is set to `$VASP_PARALLEL_CMD` by default. This means that Custodian will run VASP based on however the `VASP_PARALLEL_CMD` environment variable is set. For instance, you might include something like the following line in your submission script: `export VASP_PARALLEL_CMD="srun -N 2 --ntasks-per-node 64 -c 4 --cpu_bind=cores"`.
