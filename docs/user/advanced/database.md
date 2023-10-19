# Using a Database

Oftentimes, it is beneficial to store the results in a database for easy querying (like the example below). This is quite simple to do in quacc regardless of the workflow manager you are using by taking advantage of the numerous data store options in [maggma](https://github.com/materialsproject/maggma). For details on how to set up a Mongo database, refer to the corresponding [Maggma tutorial](https://materialsproject.github.io/maggma/getting_started/mongodb/)

![Mongo example](../../images/user/schema.gif)

=== "General Purpose"

    **Automated Approach**

    For a given recipe, you can have quacc automatically store the final output summaries in your desired database by defining a [Maggma data store](https://materialsproject.github.io/maggma/reference/stores/) in the `PRIMARY_STORE` [quacc setting](../settings/settings.md).

    For instance, let's pretend you have decided to make a [`MongoStore`](https://materialsproject.github.io/maggma/reference/stores/#maggma.stores.mongolike.MongoStore) be your database of choice. After defining or loading your Maggma store, you would call `.to_json()` to get a dictionary representation. You can then store this JSON, formatted as a string, in the `PRIMARY_STORE` global quacc setting.

    ```python
    from maggma.stores import MongoStore

    store = MongoStore(
        "my_db_name",
        "my_collection_name",
        username="my_username",
        password="my_password",
        host="localhost",
        port=27017,
    )
    print(store.to_json())  # This is the JSON string you would store in PRIMARY_STORE
    ```

    ```yaml title="~/.quacc.yaml"
    PRIMARY_STORE: '{"@module": "maggma.stores.mongolike", "@class": "MongoStore", "@version": "0.51.19", "database": "my_db_name", "collection_name": "my_collection_name", "host": "localhost", "port": 27017, "username": "my_username", "password": "my_password", "ssh_tunnel": null, "safe_update": false, "auth_source": "my_db_name", "mongoclient_kwargs": {}, "default_sort": null}'
    ```

    **Manual Approach**

    If you would prefer to store results in your database manually (perhaps because you are limited in terms of how much data you can store), you can use the [quacc.wflow.db.results_to_db][] function, as shown in the example below.

    ```python
    from maggma.stores import MongoStore
    from quacc.wflow.db import results_to_db

    # Let `results` be an output (or list of outputs) from quacc recipes

    # Define your database details
    store = MongoStore(
        "my_db_name",
        "my_collection_name",
        username="my_username",
        password="my_password",
        host="localhost",
        port=27017,
    )

    # Store the results
    results_to_db(store, results)
    ```

=== "Covalent"

    Covalent automatically stores all the inputs and outputs in an SQLite database, which you can find at the `"db_path"` when you run `covalent config`, and the results can be queried using the `#!Python ct.get_result(<dispatch ID>)` syntax. However, if you want to store the results in a different database of your choosing, you can do so quite easily.

    An example is shown below for storing the results in your custom database via the [quacc.wflow.db.covalent_to_db][] function.

    ```python
    from maggma.stores import MongoStore
    from quacc.wflow.db import covalent_to_db

    # Define your database credentials
    store = MongoStore(
        "my_db_name",
        "my_collection_name",
        username="my_username",
        password="my_password",
        host="localhost",
        port=27017,
    )

    # Store the results
    covalent_to_db(store)
    ```

=== "Jobflow"

    If you are using Jobflow to construct your workflows, it will automatically store the results in the database you defined during the [setup process](../../install/wflow_engines.md). No additional steps are needed.
