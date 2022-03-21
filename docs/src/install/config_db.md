# Database Configuration

## MongoDB Setup

### Introduction

One of the main advantages of using QuAcc is the ability to have all of your calculation inputs and outputs automatically deposited in a database of your choosing for easy querying. While this is not a requirement for using QuAcc, it is highly recommended. It is also the only way to use QuAcc with the [FireWorks](https://materialsproject.github.io/fireworks/) job management system.

### Setting up a Mongo Database

*Note*: If you already have a Mongo database, you can skip this section. Simply have your login credentials ready.

QuAcc relies on reading and writing to a Mongo Database (MongoDB). If you have never made a Mongo Database, there are many options to choose from. For new users, the easiest route is to use a cloud storage solution called [MongoDB Atlas](https://www.mongodb.com/atlas), which has a free tier. To set up your MongoDB with MongoDB Atlas, follow the instructions below.

1. Sign up for a free account on [MongoDB Atlas](https://www.mongodb.com/atlas).
2. Once logged in, select the "Build a Database" option under the "Database Deployments" section and choose the "Shared" free option.
3. Can click "Create Cluster" and enter your desired login credentials that will use to access your database. After waiting a minute or two, your cluster will be created, which is essentially a mini computer in the cloud.\
4. Go to the "Collections" tab of your cluster and create a new database. Give the database a unique name (e.g. "LastName_db") and create a collection where your QuAcc data will be stored (e.g. "quacc").
5. Retrieve your MongoDB URI, which is the address of your MongoDB cluster. You can find your database's URI by clicking the "Database" section in the sidebar and then selecting "Connect > Connect Your Application > Driver > Python > 3.11 or later" and copying the link that appears. It will be of the form `mongodb+srv://<username>:<password>@<host>/<database_name>`. Don't forget to include the <database_name> at the end, which you selected in Step 4.

## Jobflow DB Setup

To tell QuAcc where to store your calculations, you will need to make a `jobflow.yaml` file. This file will generally be formatted like the example below. Fill in the above fields with the appropriate values for your MongoDB cluster.

```yaml
JOB_STORE:
    docs_store:
      type: MongoStore
      host: <host name>
      port: 27017
      username: <username>
      password: <port>
      database: <database name>
      collection_name: <collection name>
```

If you are using a URI (e.g. as in the MongoDB Atlas instructions above), then you will instead have a `jobflow.yaml` file that instead looks like the example below. Here, you will put the full URI in the `host` field. The `username` and `password` are part of the URI and so should not be included elsewhere in the YAML file.

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
query: '{}'
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
