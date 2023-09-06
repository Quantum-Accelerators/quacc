# Settings Management

!!! Important

    The quacc configuration parameters (e.g. YAML file, environment variables) are only active on the machine where they are specified. In other words, each remote machine has its own custom set of configuration parameters.

## Modifying Quacc Settings

The default global quacc settings can be found in the [`settings.py`](https://github.com/quantum-accelerators/quacc/blob/main/src/quacc/settings.py) file. If you wish to modify any of the global quacc settings, there are several ways to do so.

### Using a YAML File

!!! Tip

    The quacc command-line interface (CLI) lets you update the YAML file from the terminal. You can run `quacc set PARAMETER value` to set a given parameter in the YAML file to the specified value. Similarly, you can run `quacc unset PARAMETER` to remove a parameter from the YAML file.

If you are planning to modify a given set of parameters for all of your calculations, the easiest way is to create a YAML file with custom settings. By default, quacc looks for this YAML file at `~/.quacc.yaml`. If you wish to store the YAML file somewhere else or with a different name, you can define the environment variable `QUACC_CONFIG_FILE` and point it to the YAML path of your choosing.

```yaml title="~/.quacc.yaml"
SCRATCH_DIR: /tmp # (1)!
CREATE_UNIQUE_WORKDIR: true # (2)!
```

1. This would set the quacc scratch directory to `/tmp` on whatever machine the calculations are run on.

2. This would ensure that each job in your quacc workflow is run in a unique, isolated working directory. This is often useful when running workflows in parallel, although some workflow engines like Covalent have their own mechanisms for this.

!!! Tip

    To set a value to the Python equivalent of `None` in a YAML file, use `null` as the value.

### Using Environment Variables

If you want to define quacc settings without writing them to a YAML file, you can instead modify the desired settings by defining individual environment variables with `QUACC` as the prefix. For instance, to modify the `SCRATCH_DIR` setting to be `/tmp`, simply define `QUACC_SCRATCH_DIR=/tmp` as a new environment variable. This approach is ideal when you want to modify the quacc settings for a subset of jobs, as the environment variable can be included in the job's submission script.

### Modifying the Global Settings in a Script

If you want to define quacc settings on-the-fly without writing them to a YAML file or using environment variables, you can do so within your script by modifying the `SETTINGS` object. This approach is ideal when you want to dynamically change a setting within a given Python process or if you're debugging in a Jupyter Notebook.

```python
from quacc import SETTINGS

SETTINGS.CREATE_UNIQUE_WORKDIR = True
```

## List of Quacc Settings

```title="quacc/settings.py"
--8<-- "../quacc/src/quacc/settings.py:settings"
```
