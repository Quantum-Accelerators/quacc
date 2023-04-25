# Database Configuration

## Introduction

One of the main advantages of using QuAcc is the ability to have all of your calculation inputs and outputs automatically deposited in a database of your choosing for easy querying. While this is not a requirement for using QuAcc, it is highly recommended.

## MongoDB Setup

*Note*: If you already have a Mongo database, you can skip this section. Simply have your login credentials ready.

QuAcc relies on reading and writing to a Mongo Database (MongoDB). If you have never made a Mongo Database, there are many options to choose from. For new users, the easiest route is to use a cloud storage solution called [MongoDB Atlas](https://www.mongodb.com/atlas), which has a free tier. To set up your MongoDB with MongoDB Atlas, follow the instructions below.

1. Sign up for a free account on [MongoDB Atlas](https://www.mongodb.com/atlas).
2. Once logged in, select the "Build a Database" option under the "Database Deployments" section and choose the "Shared" free option.
3. Select "Create Cluster" and enter your desired login credentials that will use to access your database. After waiting a minute or two, your cluster will be created, which is essentially a mini computer in the cloud.\
4. Go to the "Collections" tab of your cluster and create a new database. Give the database a unique name (e.g. "LastName_db") and create a collection where your QuAcc data will be stored (e.g. "quacc").
5. Retrieve your MongoDB URI, which is the address of your MongoDB cluster. You can find your database's URI by clicking the "Database" section in the sidebar and then selecting "Connect > Connect Your Application > Driver > Python > 3.11 or later" and copying the link that appears. It will be of the form `mongodb+srv://<username>:<password>@<host>/<database_name>`. Don't forget to include the <database_name> at the end, which you selected in Step 4.
