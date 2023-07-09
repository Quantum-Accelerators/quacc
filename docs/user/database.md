# Using a Database

Oftentimes, it is beneficial to store the results in a database for easy querying (like the example below). This is quite simple to do in quacc regardless of the workflow manager you are using by taking advantage of the numerous data store options in [maggma](https://github.com/materialsproject/maggma). Two common examples are shown here.

![Mongo example](../images/user/schema.gif)

=== "MontyDB"

    The easiest option for getting started is to use the on-disk [`MontyStore`](https://materialsproject.github.io/maggma/reference/stores/#maggma.stores.mongolike.MontyStore) option in Maggma, which acts as an interface to [MontyDB](https://github.com/davidlatwe/montydb). The only setup required is to run `pip install montydb` beforehand.

    === "General Purpose"

        For a given recipe, you can store the final output summary in your local MontyDB database using the [`quacc.util.db.results_to_db`](https://quantum-accelerators.github.io/quacc/reference/quacc/util/db.html#quacc.util.db.results_to_db) function, as shown in the example below.

        ```python
        from maggma.stores import MontyStore
        from quacc.util.db import results_to_db

        # Let `results` be an output (or list of outputs) from quacc recipes

        # Define your database credentials
        store = MontyStore(
            "quacc_results",
            database_path=".",
            database_name="quacc_db"
        )

        # Store the results
        results_to_db(store, results)
        ```

    === "Covalent"

        Covalent automatically stores all the inputs and outputs in an SQLite database, which you can find at the `"db_path"` when you run `covalent config`, and the results can be queried using the `ct.get_result(<dispatch ID>)` syntax. However, if you want to store the results in a different database of your choosing, you can do so quite easily.

        An example is shown below for storing the results in MontyDB via the [`quacc.util.db.covalent_to_db`](https://quantum-accelerators.github.io/quacc/reference/quacc/util/db.html#quacc.util.db.covalent_to_db) function.

        ```python
        from maggma.stores import MontyStore
        from quacc.util.db import covavlent_to_db

        # Define your database credentials
        store = MontyStore(
            "quacc_results",
            database_path=".",
            database_name="quacc_db"
        )

        # Store the results
        covalent_to_db(store)
        ```

    === "Jobflow"

        If you are using Jobflow to construct your workflows, it will automatically store the results in the database you defined during the [setup process](../install/wflow_engines.md).

=== "MongoDB"

    MongoDB is often preferred if you have experience with databases. For assistance with setting up a MongoDB instance of your own, refer to the ["Database Setup"](../install/config_db.md) section of the installation instructions.

    === "General Purpose"

        For a given recipe, you can store the final output summary in your database using the [`quacc.util.db.results_to_db`](https://quantum-accelerators.github.io/quacc/reference/quacc/util/db.html#quacc.util.db.results_to_db) function, as shown in the example below.

        ```python
        from maggma.stores import MongoStore
        from quacc.util.db import results_to_db

        # Let `results` be an output (or list of outputs) from quacc recipes

        # Define your database credentials
        store = MongoStore(
            "my_db_name",
            "my_collection_name",
            username="my_username",
            password="my_password",
            host="localhost",
            port=27017
        )

        # Store the results
        results_to_db(store, results)
        ```

    === "Covalent"

        Covalent automatically stores all the inputs and outputs in an SQLite database, which you can find at the `"db_path"` when you run `covalent config`, and the results can be queried using the `ct.get_result(<dispatch ID>)` syntax. However, if you want to store the results in a different database of your choosing, you can use [maggma](https://github.com/materialsproject/maggma) to do so quite easily.

        An example is shown below for storing the results in a MongoDB via the [`quacc.util.db.covalent_to_db`](https://quantum-accelerators.github.io/quacc/reference/quacc/util/db.html#quacc.util.db.covalent_to_db) function.

        ```python
        from maggma.stores import MongoStore
        from quacc.util.db import covavlent_to_db

        # Define your database credentials
        store = MongoStore(
            "my_db_name",
            "my_collection_name",
            username="my_username",
            password="my_password",
            host="localhost",
            port=27017
        )

        # Store the results
        covalent_to_db(store)
        ```

    === "Jobflow"

        If you are using Jobflow to construct your workflows, it will automatically store the results in the database you defined during the [setup process](../install/wflow_engines.md).
