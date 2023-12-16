# Settings Management

!!! Important

    The quacc configuration parameters (e.g. YAML file, environment variables) are only active on the machine where they are specified. In other words, each remote machine has its own custom set of configuration parameters.

## Modifying Quacc Settings

The default global quacc settings can be found in the [`settings.py`](https://github.com/quantum-accelerators/quacc/blob/main/src/quacc/settings.py) file. If you wish to modify any of the global quacc settings, there are several ways to do so.

### Using a YAML File

If you are planning to modify a given set of parameters for all of your calculations, the easiest way is to create a YAML file with custom settings. By default, quacc looks for this YAML file at `~/.quacc.yaml`. If you wish to store the YAML file somewhere else or with a different name, you can define the environment variable `QUACC_CONFIG_FILE` and point it to the YAML path of your choosing.

```yaml title="~/.quacc.yaml"
SCRATCH_DIR: $SCRATCH # (1)!
CREATE_UNIQUE_DIR: true # (2)!
WORKFLOW_ENGINE: # (3)!
```

1. In YAML, you do not need to put quotation marks around strings. In this case, the quacc scratch directory would be set to the environment variable `$SCRATCH` on whatever machine the calculations are run on.

2. In YAML, booleans are lowercase. In this case, the setting would ensure that each job in your quacc workflow is stored in a unique, isolated directory.

3. In YAML, a blank value is interpreted as `None` in Python. A value of `null` is equivalent to `None` too. In this case, no workflow engine would be used.

!!! Tip "Command-Line Interface"

    The quacc command-line interface (CLI) lets you update the YAML file from the terminal. You can run `quacc set PARAMETER value` to set a given parameter in the YAML file to the specified value. Similarly, you can run `quacc unset PARAMETER` to remove a parameter from the YAML file.

### Using Environment Variables

If you want to define quacc settings without writing them to a YAML file, you can instead modify the desired settings by defining individual environment variables with `QUACC` as the prefix. For instance, to modify the `SCRATCH_DIR` setting to be `$SCRATCH`, simply define `QUACC_SCRATCH_DIR=$SCRATCH` as a new environment variable. This approach is ideal when you want to modify the quacc settings for a subset of jobs, as the environment variable can be included in the job's submission script.
