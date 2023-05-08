# Fireworks Configuration

## Introduction

Follow the instructions to configure Quacc to work with Jobflow+FireWorks.

## Jobflow DB Setup

To tell QuAcc where to store your calculations, you will need to make a `jobflow.yaml` file. This file will generally be formatted like the example below. Fill in the above fields with the appropriate values for your MongoDB cluster.

```yaml
JOB_STORE:
    docs_store:
      type: MongoStore
      host: <host name>
      port: 27017
      username: <username>
      password: <password>
      database: <database name>
      collection_name: <collection name>
```

If you are using a URI (e.g. as in the MongoDB Atlas instructions above), then you will instead have a `jobflow.yaml` file that looks like the example below. Here, you will put the full URI in the `host` field. The `username` and `password` are part of the URI and so should not be included elsewhere in the YAML file.

```yaml
JOB_STORE:
    docs_store:
      type: MongoStore
      host: <URI>
      port: 27017
      database: <database name>
      collection_name: <collection name>
```

Finally, you will need to define a `JOBFLOW_CONFIG_FILE` environment variable pointing to the file you made. For instance, in your `~/.bashrc` file, add the following line:
`export JOBFLOW_CONFIG_FILE="/path/to/my/jobflow.yaml"`.

When a QuAcc calculation completes, all the data will be stored in the database you have specified above. If needed, additional information about setting up your `jobflow.yaml` file can be found in the official [jobflow documentation](https://materialsproject.github.io/jobflow/jobflow.settings.html).

## FireWorks DB Setup

If you plan to use FireWorks, you will also need to make a few configuration files: `FW_config.yaml`, `my_fworker.yaml`, `my_launchpad.yaml`, and `my_qadapter.yaml`.

To begin, make a directory called `fw_config` where you will store the above four files.

### FW Config File

For the `FW_config.yaml`, you can use the following template. Make sure to update the path to the `fw_config` folder you made earlier.

```yaml
CONFIG_FILE_DIR: </path/to/fw_config>
QUEUE_UPDATE_INTERVAL: 2
```

### FWorker

For the `my_fworker.yaml`, you can use the following template. You do not need to make any modifications.

```yaml
name: quacc_fworker
category: ''
query: {}
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

Assuming you plan to use a queuing system for your compute jobs, you will need to make a `my_qadapter.yaml` file. For this, you will need to follow the instructions in the [FireWorks documentation](https://materialsproject.github.io/fireworks/qadapter_programming.html) for your specific job scheduling system. An example `my_qadapter.yaml` file is shown below for SLURM.

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
            export VASP_PARALLEL_CMD="srun -N 2 --ntasks-per-node=24"
```

In the above example, you would need to change the path in the `rocket_launch` field to the correct path to your `my_fworker.yaml`. The nodes, walltime, account, and qos are the corresponding parameters for your queuing system. Finally, anything in the `pre_rocket` field will be executed before the job begins running. It is a good place to load modules and set environment variables. A representative example has been provided above.

Finally, you will need to define a `FW_CONFIG_FILE` environment variable pointing to the `FW_config.yaml` file you made. For instance, in your `~/.bashrc` file, add the following line:
`export FW_CONFIG_FILE="/path/to/config/fw_config/FW_config.yaml"`.


### Initialization

To check that everything is working right with FireWorks, run lpad reset if you haven’t run it before to ensure there is a connection to the database.