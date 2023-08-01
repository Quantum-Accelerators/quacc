# Quacc Settings

## Modifying Quacc Settings

By default, the global quacc settings can be found in the [`settings.py`](https://github.com/quantum-accelerators/quacc/blob/main/src/quacc/settings.py) file. If you wish to modify any of the global quacc settings, there are two ways to do so.

### Using a YAML File

The easiest way is to create a custom YAML file with custom settings. For instance, it might look something like the following if you wanted to modify the `SCRATCH_DIR`, `GZIP_FILES`, and `CREATE_UNIQUE_WORKDIR` settings.

```yaml title="quacc.yaml"
SCRATCH_DIR: /tmp # (1)!
GZIP_FILES: false # (2)!
CREATE_UNIQUE_WORKDIR: true # (3)!
```

1. This would set the quacc scratch directory to `/tmp` on whatever machine the calculations are run on.

2. This would disable the automatic gzip-ing of all files stored in the quacc results directory.

3. This would ensure that each job in your quacc workflow is run in a unique, isolated working directory. This is often useful when running workflows in parallel, although some workflow engines like Covalent have their own mechanisms for this.

By default, quacc looks for this YAML file at `~/quacc.yaml`. If you wish to store the `quacc.yaml` file somewhere else, you can define the environment variable `QUACC_CONFIG_FILE` and point it to the `quacc.yaml` path of your choosing.

### Using Environment Variables

If you want to define quacc settings on-the-fly without writing them to a YAML file, you can instead modify the desired settings by defining individual environment variables with `QUACC` as the prefix. For instance, to modify the `SCRATCH_DIR` setting to be `/tmp`, simply define `QUACC_SCRATCH_DIR=/tmp` as a new environment variable.
