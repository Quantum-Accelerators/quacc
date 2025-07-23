# Modifying Settings

!!! Important

    The quacc configuration parameters (e.g. YAML file, environment variables) are only active on the machine where they are specified. In other words, each remote machine has its own custom set of configuration parameters.

The default global quacc settings can be found in the [`settings.py`](https://github.com/quantum-accelerators/quacc/blob/main/src/quacc/settings.py) file. If you wish to modify any of the global quacc settings, there are several ways to do so.

## Using a YAML File

!!! Note "Command-Line Interface"

    The quacc command-line interface (CLI) lets you update the YAML file from the terminal. You can run `quacc set PARAMETER value` to set a given parameter in the YAML file to the specified value. Similarly, you can run `quacc unset PARAMETER` to remove a parameter from the YAML file.

If you are planning to modify a given set of parameters for all of your calculations, the easiest way is to create a YAML file with custom settings. By default, quacc looks for this YAML file at `~/.quacc.yaml`. If you wish to store the YAML file somewhere else or with a different name, you can define the environment variable `QUACC_CONFIG_FILE` and point it to the YAML path of your choosing.

```yaml title="~/.quacc.yaml"
SCRATCH_DIR: /path/to/my/scratch/dir # (1)!
CREATE_UNIQUE_DIR: false # (2)!
WORKFLOW_ENGINE: None # (3)!
```

1. In YAML, you do not need to put quotation marks around strings. The YAML file also supports environment variable-based paths, like `$SCRATCH`.

2. In YAML, booleans are lowercase.

3. In YAML, a blank value or `null` is interpreted as `None` in Python. We also support "None" for convenience.

!!! Tip "When is This Method Ideal?"

    This approach is ideal when you want to change a setting that applies to most or all of your calculations. If you wish to have multiple sets of YAML files, you can dynamically change which file to use by setting the `QUACC_CONFIG_FILE` environment variable accordingly.

## Using Environment Variables

If you want to define quacc settings without writing them to a YAML file, you can instead modify the desired settings by defining individual environment variables with `QUACC` as the prefix. The environment variable takes precedence over any value specified in the YAML file. Most simple field types (e.g. `int`, `bool`, `float`, `str`) will be automatically inferred from the environment variable. For more complex types, such as lists or dictionaries, refer to the corresponding section in the [pydantic-settings](https://docs.pydantic.dev/latest/concepts/pydantic_settings/#parsing-environment-variable-values) documentation. We use `__` as the delimiter for nested settings.

To achieve the same results as the aforementioned YAML file, you would define the following environment variables:

```bash
export QUACC_SCRATCH_DIR=/path/to/my/scratch/dir
export QUACC_CREATE_UNIQUE_DIR=False
export QUACC_WORKFLOW_ENGINE=None
```

!!! Tip "When is This Method Ideal?"

    This approach is ideal when you want to dynamically modify the quacc settings for a given set of calculations, as the environment variable can be modified on-the-fly (e.g. in the job's submission script) without modifying the YAML file that is read by all other calculations.

!!! Warning "Be Careful!"

    If you try to modify environment variables within Python using `os.environ`, make sure you do this before you make any `quacc`-related imports.

## Modifying the Global Settings in a Script

If you want to define quacc settings on-the-fly without writing them to a YAML file or using environment variables, you can do so using the context handler function [quacc.settings.change_settings][] as follows:

```python
from ase.build import bulk
from quacc import change_settings
from quacc.recipes.emt.core import relax_job

atoms = bulk("Cu")

with change_settings({"GZIP_FILES": False}):
    result = relax_job(atoms)
```

??? Warning "Using a Workflow Engine?"

    When deploying calculations via a workflow engine, changes to in-memory global variables on the local machine will not be reflected on the remote machine. While the YAML file and environment variable approaches are preferred, if you wish to make changes from within your script then this should be done via a custom `settings_swap` keyword argument that is supported by the `@job` decorator.

    ```python
    from quacc import job


    @job(settings_swap={"GZIP_FILES": False})  # (1)!
    def add(a, b):
        return a + b
    ```

    1. This is the same as doing

        ```python
        from quacc import change_settings, job


        @job
        def add(a, b):
            with change_settings({"GZIP_FILES": False}):
                return a + b
        ```

    If using a pre-made `@job`, you can simply redecorate it so that it supports your custom settings:

    ```python
    from ase.build import bulk
    from quacc import redecorate, job
    from quacc.recipes.emt.core import relax_job

    atoms = bulk("Cu")
    relax_job_ = redecorate(relax_job, job(settings_swap={"GZIP_FILES": False}))
    results = relax_job_(atoms)
    ```

!!! Tip "When is This Method Ideal?"

    This approach is ideal for fine-tuned modifications to settings within your workflow.
