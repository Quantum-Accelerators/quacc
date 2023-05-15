# Database Usage

Oftentimes, it is beneficial to store the results in a database for easy querying. There are many ways to achieve this with QuAcc depending on your preferences.

## With Covalent as the Workflow Manager

Covalent automatically stores all the inputs and outputs in an SQLite database, which you can find at the `"db_path"` when you run `covalent config`, and the results can be queried using the `ct.get_result(<dispath ID>)` syntax. However, if you want to store the results in a different database of your choosing, you can use [`maggma`](https://github.com/materialsproject/maggma) to do so quite easily. An example is shown below for storing the results in a MongoDB. For assistance with setting up a MongoDB of your own, refer to the "Optional: MongoDB Setup" section of the installation instructions.

```python
import covalent as ct
from maggma.stores import MongoStore

# Connect to the database
database = "my_db"
collection_name = "my_collection"
store = MongoStore(database, collection_name, host="localhost", port=27017, username="my_username", password="my_password")
store.connect()

# Fetch the results
results_dir = ct.get_config()["dispatcher"]["results_dir"]
docs = []
for dispatch_id in os.listdir(results_dir):
    result = ct.get_result(dispatch_id).result
    docs.append({"dispatch_id":dispatch_id, "result": result})

# Store the results
with store:
    store.update(docs, key="dispatch_id")

# Close the database connection
store.close()
```

## With Jobflow as the Workflow Manager

If you are using Jobflow to construct your workflows, it will automatically store the results in the database you defined during the setup process.