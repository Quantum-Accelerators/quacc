# Database Options

## Why Use a Database?

In most cases, it is beneficial to store the results in a database for easy querying (like the example below). This is quite simple to do in quacc regardless of the workflow manager you are using by taking advantage of the numerous data store options in [maggma](https://github.com/materialsproject/maggma).

![Mongo example](../../images/user/schema.gif)

## Choosing a Data Store Option

Some workflow engines supported by quacc, such as Covalent and Jobflow, provide a database for storing results by default. In these cases, you can simply use the database provided by the workflow engine. However, if you are using a workflow engine that does not provide a database or not using a workflow engine at all, you may wish to use one directly supported by quacc.

The supported data store options are all those [provided by maggma](https://materialsproject.github.io/maggma/getting_started/stores/#list-of-stores). At the time of writing, two of the most popular options are:

- A [`MongoStore`](https://materialsproject.github.io/maggma/reference/stores/#maggma.stores.mongolike.MemoryStore) for interfacing with MongoDB. This is generally recommended, especially if your compute nodes have access to a network.
- A [`MontyStore`](https://materialsproject.github.io/maggma/reference/stores/#maggma.stores.mongolike.MontyStore) for interfacing with a filesystem-based storage option that mimics MongoDB. This is generally recommended if setting up or connecting to a MongoDB instance is not practical.
