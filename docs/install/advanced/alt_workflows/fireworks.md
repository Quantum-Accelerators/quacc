# FireWorks Setup

## Introduction

[FireWorks](https://github.com/materialsproject/fireworks) can be used to dispatch and monitor workflows, particularly those made using Jobflow. Follow the instructions below for how to use FireWorks with quacc. For additional details, refer to the full [Fireworks documentation](https://github.com/materialsproject/fireworks).

## Installation

To install quacc with support for FireWorks, run `pip install fireworks`.

## FireWorks DB Setup

If you plan to use FireWorks to dispatch your Jobflow workflows, you will also need to make a few configuration files: `FW_config.yaml`, `my_fworker.yaml`, `my_launchpad.yaml`, and `my_qadapter.yaml`. To begin, make a directory called `fw_config` where you will store the four files described in greater detail below. The directory structure will look like the following:

```text
fw_config
├── FW_config.yaml
├── my_fworker.yaml
├── my_launchpad.yaml
└── my_qadapter.yaml
```

### FW Config File

For the `FW_config.yaml`, you can use the following template. Make sure to update the path to the `fw_config` folder where the file resides.

```yaml
CONFIG_FILE_DIR: </path/to/fw_config>
QUEUE_UPDATE_INTERVAL: 2
```

### FWorker

For the `my_fworker.yaml`, you can use the following template. You do not need to make any modifications.

```yaml
name: quacc_fworker
category: ""
query: "{}"
```

### Launchpad

For the `my_launchpad.yaml`, you can use the following template. Replace the entries in `<>` with the appropriate values for your Mongo database.

```yaml
host: <host name>
port: 27017
name: <database name>
username: <username>
password: <password>
logdir: null
Istrm_lvl: DEBUG
user_indices: []
wf_user_indices: []
```

If you are accessing your MongoDB via a URI (e.g. as with MongoDB Atlas), then you will use the following `my_launchpad.yaml` template instead.

```yaml
host: <URI>
port: 27017
name: <database name>
uri_store: true
logdir: null
Istrm_lvl: DEBUG
user_indices: []
wf_user_indices: []
```

### QAdapter

Assuming you plan to use a queuing system for your compute jobs, you will need to make a `my_qadapter.yaml` file. For this, you will need to follow the instructions in the [FireWorks documentation](https://materialsproject.github.io/fireworks/qadapter_programming.html) for your specific job scheduling system. An example `my_qadapter.yaml` file is shown below for Slurm.

```yaml
_fw_name: CommonAdapter
_fw_q_type: SLURM
rocket_launch: rlaunch -w /path/to/fw_config/my_fworker.yaml singleshot
nodes: 2
walltime: 00:30:00
account: <account>
job_name: quacc_firework
qos: regular
pre_rocket: |
  module load vasp
  export QUACC_VASP_PARALLEL_CMD="srun -N 2 --ntasks-per-node=24 --cpu_bind=cores"
```

In the above example, you would need to change the path in the `rocket_launch` field to the correct path to your `my_fworker.yaml`. The nodes, walltime, account, and qos are the corresponding parameters for your queuing system. Finally, anything in the `pre_rocket` field will be executed before the job begins running. It is a good place to load modules and set environment variables. A representative example has been provided above.

Finally, you will need to define a `FW_CONFIG_FILE` environment variable pointing to the `FW_config.yaml` file you made. For instance, in your `~/.bashrc` file, add the following line:
`export FW_CONFIG_FILE="/path/to/config/fw_config/FW_config.yaml"`.

### Database Initialization

```{warning}
Running `lpad reset` will clear your FireWorks launchpad, so only use this command if you are a new user.
```

To check that everything is working right with FireWorks, run `lpad reset` to ensure there is a connection to the database.
