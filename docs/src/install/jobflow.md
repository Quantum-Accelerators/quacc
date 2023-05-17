# Optional: Jobflow Setup

## Introduction

While we recommend using [Covalent](https://github.com/AgnostiqHQ/covalent) as the workflow manager with Quacc, it is not the only option for you to consider.

If you would prefer to use [Jobflow](https://github.com/materialsproject/jobflow) to write your workflows and/or [FireWorks](https://github.com/materialsproject/fireworks) to manage them, follow the instructions below. For additional details, refer to the full [Jobflow documentation](https://materialsproject.github.io/jobflow/) and [FireWorks documentation](https://materialsproject.github.io/fireworks/).

## MongoDB Setup

Jobflow and FireWorks both require the use of a MongoDB database to store calculation results. If you haven't done so already, first create a Mongo database as described in the ["Optional: MongoDB Configuration"](covalent.md) section.

## Jobflow DB Setup

If you plan to use Jobflow to write your workflows, you will need to make a `jobflow.yaml` file. This file will generally be formatted like the example below. Fill in the above fields with the appropriate values for your MongoDB cluster.

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

If you are using a URI (e.g. as is common with MongoDB Atlas), then you will instead have a `jobflow.yaml` file that looks like the example below. Here, you will put the full URI in the `host` field. The `username` and `password` are part of the URI and so should not be included elsewhere in the YAML file.

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

When a Quacc calculation completes, all the data will be stored in the database you have specified above. If needed, additional information about setting up your `jobflow.yaml` file can be found in the official [jobflow documentation](https://materialsproject.github.io/jobflow/jobflow.settings.html).
