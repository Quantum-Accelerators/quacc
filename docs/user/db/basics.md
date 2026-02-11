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

    Then, to store results in your database, you can use the following script:

    ```python
    from maggma.stores import MongoStore
    from monty.json import jsanitize

    # Let `results` be of type `list[dict]` containing outputs from quacc recipes

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
    sanitized_results = [
        jsanitize(result, enum_values=True, recursive_msonable=True) for result in results
    ]

    for result in sanitized_results:
        result["uuid"] = str(uuid.uuid4())

    with store:
        store.update(sanitized_results, key="uuid")
    ```

=== "Prefect"

    To use a database with Prefect, read the [Database](https://docs.prefect.io/latest/concepts/database/) section of the Prefect documentation as well as how to [persist results](https://docs.prefect.io/latest/concepts/results/#persisting-results).

=== "Jobflow"

    If you are using Jobflow to construct your workflows, it will automatically store the results in the database you defined during the [setup process](../../install/wflow_engines.md). No additional steps are needed.
