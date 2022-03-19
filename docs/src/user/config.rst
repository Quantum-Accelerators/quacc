
=============
Configuration
=============

Database Setup
==============
One of the main advantages of using QuAcc is the ability to have all of your calculation inputs and outputs automatically deposited in a database of your choosing for easy querying. While this is not a requirement for using QuAcc, it is highly recommended. It is also the only way to use QuAcc with the [FireWorks](https://materialsproject.github.io/fireworks/) job management system.

QuAcc relies on reading and writing to a Mongo Database (MongoDB). If you have never made a Mongo Database, there are many options to choose from. If you are new to MongoDB, the easiest route is to use a cloud storage solution called [MongoDB Atlas](https://www.mongodb.com/atlas), which has a free tier. To set up your MongoDB with MongoDB Atlas, sign up for a free account. Once logged in, select the "Build a Database" option under the "Database Deployments" section and choose the "Shared" free option. Then can click "Create Cluster" and enter your desired login credentials that will use to access your database. After waiting a minute or two, your cluster will be created, which is essentially a mini computer in the cloud. The last step is to go to the "Collections" tab of your cluster and create a new database. Give the database a unique name (e.g. "LastName_db") and create a collection where your QuAcc data will be stored (e.g. "quacc"). If using MongoDB Atlas, the last step is to find your MongoDB URI, which is the address of your MongoDB cluster. You can find your database's URI by clicking the "Database" section in the sidebar and then selecting "Connect > Connect Your Application > Driver > Python > 3.11 or later" and copying the link that appears. It will be of the form "mongodb+srv://<username>:<password>@<host>/<database_name>"

If you already are using MongoDB and have a database available, you do not need to follow the above steps. Simply have your login credentials ready.

Jobflow DB Setup
================
If you completed the above database setup, then QuAcc can store all the input and output data for each calculation in the database you made. While this is optional, it is the best way to take advantage of the features in QuAcc.

To tell QuAcc where to store your calculations, you will need to make a `jobflow.yaml` file. This file will be formatted as follows:
```yaml
JOB_STORE:
    docs_store:
      type: MongoStore
      host: <host name>
      port: <port>
      username: <username>
      password: <port>
      database: <database name>
      collection_name: <collection name>
```
Fill in the above fields with the appropriate values for your MongoDB cluster.

If you are using MongoDB Atlas or any other cloud service with a provided URI, then you will instead have a `jobflow.yaml` file that looks like
```yaml
JOB_STORE:
    docs_store:
      type: MongoStore
      host: <URI>
      port: <port>
      database: <database name>
      collection_name: <collection name>
```
Here, you will put the full URI in the `host` field. The username and password are part of the URI and so should not be included elsewhere in the YAML file.

Additional information about setting up your `jobflow.yaml` file can be found in the official [jobflow documentation](https://materialsproject.github.io/jobflow/jobflow.settings.html).

With your `jobflow.yaml` file completed, you will need to define a `JOBFLOW_CONFIG_FILE` environment variable pointing to the file you made. For instance, in your `~/.bashrc` file, add the following line:
`export JOBFLOW_CONFIG_FILE="/path/to/my/jobflow.yaml"`.

When a QuAcc calculation completes, all the data will be stored in the database you have specified above.

FireWorks DB Setup
==================
If you plan to use FireWorks, you will also need to make a few configuration files: `FW_config.yaml`, `my_fworker.yaml`, `my_launchpad.yaml`, and `my_qadapter.yaml`.

Make a directory called `fw_config` where you will store the above four files.

For the `FW_config.yaml`, you can use the following template:
```yaml
CONFIG_FILE_DIR: /path/to/fw_config
QUEUE_UPDATE_INTERVAL: 2
```
Make sure to update the path to the `fw_config` folder you just made.

For the `my_fworker.yaml`, you can use the following template:
```yaml
name: quacc_fworker
category: ''
query: '{}'
```

For the `my_launchpad.yaml`, you can use the following template:
```yaml
host: <host name>
port: <port>
name: <database name>
username: <username>
password: <password>
logdir: null
Istrm_lvl: DEBUG
user_indices: []
wf_user_indices: []
```
Replace the entries in <> with the appropriate values for your Mongo database. If you don't know your `port` value, you can likely use 27017.

If you are accessing your MongoDB via a URI (e.g. as with MongoDB Atlas), then you will use the following `my_launchpad.yaml` template instead:
```yaml
host: <URI>
port: <port>
name: <database name>
uri_store: true
logdir: null
Istrm_lvl: DEBUG
user_indices: []
wf_user_indices: []
```

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
            export OMP_PROC_BIND=true
            export OMP_PLACES=threads
            export OMP_NUM_THREADS=1
```

In the above example, you would need to change the path in the `rocket_launch` field to the correct `fw_config` path. The nodes, walltime, account, and qos are the corresponding parameters for your queuing system. Finally, anything in the `pre_rocket` field will be executed before the job begins running. It is a good place to load modules and set environment variables. A representative example has been provided above.

Finally, you will need to define a `FW_CONFIG_FILE` environment variable pointing to the `FW_config.yaml` file you made. For instance, in your `~/.bashrc` file, add the following line:
`export FW_CONFIG_FILE="/path/to/config/fw_config/FW_config.yaml"`.
