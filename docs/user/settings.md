# Quacc Settings

## Modifying Quacc Settings

The default global quacc settings can be found in the [`settings.py`](https://github.com/quantum-accelerators/quacc/blob/main/src/quacc/settings.py) file. If you wish to modify any of the global quacc settings, there are several ways to do so.

### Using a YAML File

If you are planning to modify parameters for most of your calculations, the easiest way is to create a YAML file with custom settings. For instance, it might look something like the following if you wanted to modify the `WORKFLOW_ENGINE`, `SCRATCH_DIR`, and `CREATE_UNIQUE_WORKDIR` settings for all your workflows.

```yaml title="~/.quacc.yaml"
WORKFLOW_ENGINE: covalent # (1)!
SCRATCH_DIR: /tmp # (2)!
CREATE_UNIQUE_WORKDIR: true # (3)!
```

1. This would set the workflow engine to Covalent. If you wish to disable the use of a workflow engine, set this to `null`.

2. This would set the quacc scratch directory to `/tmp` on whatever machine the calculations are run on.

3. This would ensure that each job in your quacc workflow is run in a unique, isolated working directory. This is often useful when running workflows in parallel, although some workflow engines like Covalent have their own mechanisms for this.

By default, quacc looks for this YAML file at `~/.quacc.yaml`. If you wish to store the YAML file somewhere else or with a different name, you can define the environment variable `QUACC_CONFIG_FILE` and point it to the YAML path of your choosing.

### Using Environment Variables

If you want to define quacc settings on-the-fly without writing them to a YAML file, you can instead modify the desired settings by defining individual environment variables with `QUACC` as the prefix. For instance, to modify the `SCRATCH_DIR` setting to be `/tmp`, simply define `QUACC_SCRATCH_DIR=/tmp` as a new environment variable.

### Modifying the Global Settings in a Script

If you want to define quacc settings on-the-fly without writing them to a YAML file or using environment variables, you can do so within your script by modifying the `SETTINGS` object as follows:

```python
from quacc import SETTINGS

SETTINGS.CREATE_UNIQUE_WORKDIR = True
```

!!! Important

    This approach should work as expected for all settings except `WORKFLOW_ENGINE`, which must be defined before any importing any recipes. Therefore, it is recommended to use the YAML file or environment variable approach for specifying the workflow engine.

## List of Quacc Settings

```title="quacc/settings.py"
--8<-- "../quacc/src/quacc/settings.py:settings"
```
