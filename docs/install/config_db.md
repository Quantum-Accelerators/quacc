# Database Setup

## Introduction

Some users may wish to use quacc in way that ensures all calculation inputs and outputs are stored in an easily queryable database. While not a requirement, this can be readily achieved through the use of the [maggma](https://github.com/materialsproject/maggma) package. Maggma has many options for [data stores](https://materialsproject.github.io/maggma/reference/stores/), but we will only describe how to set up a MongoDB store since it is the most common.

## MongoDB Setup

!!! Note

    Users of NERSC HPC machines can instead [request a database](https://docs.nersc.gov/services/databases/) directly from NERSC staff.

For new users, the easiest route to create a Mongo database is to use a cloud storage solution called [MongoDB Atlas](https://www.mongodb.com/atlas), which has a free tier. To set up your MongoDB with MongoDB Atlas, follow the instructions below:

1. Sign up for a free account on [MongoDB Atlas](https://www.mongodb.com/atlas).
2. Once logged in, select the "Build a Database" option under the "Database Deployments" section and choose the "Shared" free option.
3. Select "Create Cluster" and enter your desired login credentials that will use to access your database. After waiting a minute or two, your cluster will be created, which is essentially a mini computer in the cloud.
4. Go to the "Collections" tab of your cluster and create a new database. Give the database a unique name (e.g. "LastName_db") and create a collection where your quacc data will be stored (e.g. "quacc").
5. Retrieve your MongoDB URI, which is the address of your MongoDB cluster. You can find your database's URI by clicking the "Database" section in the sidebar and then selecting "Connect > Connect Your Application > Driver > Python > 3.11 or later" and copying the link that appears. It will be of the form `mongodb+srv://<username>:<password>@<host>/<database_name>`. Don't forget to include the <database_name> at the end, which you selected in Step 4.
