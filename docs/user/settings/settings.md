# Modifying Settings

!!! Important

    The quacc configuration parameters (e.g. YAML file, environment variables) are only active on the machine where they are specified. In other words, each remote machine has its own custom set of configuration parameters.

The default global quacc settings can be found in the [`settings.py`](https://github.com/quantum-accelerators/quacc/blob/main/src/quacc/settings.py) file. If you wish to modify any of the global quacc settings, there are several ways to do so.

## Using a YAML File

!!! Note "Command-Line Interface"

    The quacc command-line interface (CLI) lets you update the YAML file from the terminal. You can run `quacc set PARAMETER value` to set a given parameter in the YAML file to the specified value. Similarly, you can run `quacc unset PARAMETER` to remove a parameter from the YAML file.

If you are planning to modify a given set of parameters for all of your calculations, the easiest way is to create a YAML file with custom settings. By default, quacc looks for this YAML file at `~/.quacc.yaml`. If you wish to store the YAML file somewhere else or with a different name, you can define the environment variable `QUACC_CONFIG_FILE` and point it to the YAML path of your choosing.

```yaml title="~/.quacc.yaml"
SCRATCH_DIR: $SCRATCH # (1)!
CREATE_UNIQUE_DIR: false # (2)!
WORKFLOW_ENGINE: None # (3)!
```

1. In YAML, you do not need to put quotation marks around strings. In this case, the quacc scratch directory would be set to the environment variable `$SCRATCH` on whatever machine the calculations are run on.

2. In YAML, booleans are lowercase.

3. In YAML, a blank value or `null` is interpreted as `None` in Python. We also support "None" for convenience.

??? Tip "When is This Method Ideal?"

    This approach is ideal when you want to change a setting that applies to most or all of your calculations.

## Using Environment Variables

If you want to define quacc settings without writing them to a YAML file, you can instead modify the desired settings by defining individual environment variables with `QUACC` as the prefix. The environment variable takes precedence over any value specified in the YAML file. Most simple field types (e.g. `int`, `bool`, `float`, `str`) will be automatically inferred from the environment variable. To achieve the same results as the aforementioned YAML file, you would define the following environment variables:

```bash
export QUACC_SCRATCH_DIR=$SCRATCH
export QUACC_CREATE_UNIQUE_DIR=False
export QUACC_WORKFLOW_ENGINE=None
```

For more complex types, such as lists or dictionaries, refer to the corresponding section in the [pydantic-settings](https://docs.pydantic.dev/latest/concepts/pydantic_settings/#parsing-environment-variable-values) documentation. We use `__` as the delimiter for nested settings.

??? Tip "When is This Method Ideal?"

    This approach is ideal when you want to dynamically modify the quacc settings for a given set of calculations, as the environment variable can be modified on-the-fly (e.g. in the job's submission script) without modifying the YAML file that is read by all other calculations.

## Modifying the Global Settings in a Script

If you want to define quacc settings on-the-fly without writing them to a YAML file or using environment variables, you can do so within your script by modifying the global `SETTINGS` object.

```python
from quacc import SETTINGS

SETTINGS.RESULTS_DIR = "/new/path/to/store/results"
```

??? Tip "When is This Method Ideal?"

    This approach is ideal for debugging scenarios, such as when using a Jupyter Notebook. Generally, it should not be used when deploying calculations via a workflow engine, as changes to in-memory global variables on the local machine will not be reflected on the remote machine.
