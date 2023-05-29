# Jobflow Setup

## Introduction

While we recommend using [Covalent](https://github.com/AgnostiqHQ/covalent) as the workflow manager with Quacc, it is not the only option for you to consider. If you would prefer to use [Jobflow](https://github.com/materialsproject/jobflow) to write your workflows and/or [FireWorks](https://github.com/materialsproject/fireworks) to manage them, follow the instructions below. For additional details, refer to the full [Jobflow documentation](https://materialsproject.github.io/jobflow/).

## MongoDB Setup

Jobflow and FireWorks both require the use of a database to store calculation results. If you haven't done so already, first create a Mongo database as described in the ["MongoDB Setup"](config_db.md) section.

## Jobflow DB Setup

If you plan to use Jobflow to write your workflows, you will need to make a `jobflow.yaml` file. This file will generally be formatted like the example below. Fill in the fields with the appropriate values for your MongoDB cluster, which is where all your calculation inputs and outputs will be stored.

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

If you are using a URI (as is common with MongoDB Atlas), then you will instead have a `jobflow.yaml` file that looks like the example below. Here, you will put the full URI in the `host` field. The `username` and `password` are part of the URI and so should not be included elsewhere in the YAML file.

```yaml
JOB_STORE:
  docs_store:
    type: MongoStore
    host: <URI>
    port: 27017
    database: <database name>
    collection_name: <collection name>
```

You will then need to define a `JOBFLOW_CONFIG_FILE` environment variable pointing to the file you made. For instance, in your `~/.bashrc` file, add the following line:
`export JOBFLOW_CONFIG_FILE="/path/to/my/jobflow.yaml"`.
