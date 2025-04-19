# Setup and Storing Results

Here, we describe how to set up quacc with a database of your choosing.

=== "General Purpose"

    You can store the output of quacc jobs and flows in your database of choice by defining a [Maggma data store](https://materialsproject.github.io/maggma/reference/stores/).

    For instance, let's pretend you have decided to make a [`MongoStore`](https://materialsproject.github.io/maggma/reference/stores/#maggma.stores.mongolike.MongoStore) be your database of choice. Instantiating that class might look like the following

    ```python
    from maggma.stores import MongoStore

    store = MongoStore(
        database="my_db_name",
        collection_name="my_collection_name",
        username="my_username",
        password="my_password",
        host="my_hostname",
        port=27017,
    )
    ```

    Then, to store results in your database, you can use the [quacc.wflow_tools.db.results_to_db][] function, as shown in the example below.

    ```python
    from maggma.stores import MongoStore
    from quacc.wflow_tools.db import results_to_db

    # Let `results` be an output (or list of outputs) from quacc recipes

    # Define your database details
    store = MongoStore(
        database="my_db_name",
        collection_name="my_collection_name",
        username="my_username",
        password="my_password",
        host="my_hostname",
        port=27017,
    )

    # Store the results
    results_to_db(store, results)
    ```

=== "Covalent"

    Covalent automatically stores all the inputs and outputs in an SQLite database, which you can find at the `"db_path"` when you run `covalent config`, and the results can be queried using the `#!Python ct.get_result(<dispatch ID>)` syntax.

=== "Prefect"

    To use a database with Prefect, read the [Database](https://docs.prefect.io/latest/concepts/database/) section of the Prefect documentation as well as how to [persist results](https://docs.prefect.io/latest/concepts/results/#persisting-results).

=== "Jobflow"

    If you are using Jobflow to construct your workflows, it will automatically store the results in the database you defined during the [setup process](../../install/wflow_engines.md). No additional steps are needed.
