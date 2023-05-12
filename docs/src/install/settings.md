# Settings

##  Modifying QuAcc Settings

By defualt, the global QuAcc settings can be found in the `settings.py` file [here](https://github.com/arosen93/quacc/blob/main/quacc/settings.py). If you wish to modify any of the global QuAcc settings, there are two ways to do so.

### Method 1

The easiest way is to create a custom YAML file with custom settings. For instance, it might look something like the following if you wanted to modify the `SCRATCH_DIR` and `GZIP_FILES` settings.

```yaml
SCRATCH_DIR: /tmp
GZIP_FILES: false
```

By default, QuAcc looks for this YAML file at `~/quacc.yaml`. If you wish to store the `quacc.yaml` file somewhere else, you can define the environment variable `QUACC_CONFIG_FILE` and point it to the `quacc.yaml` path of your choosing.

### Method 2

If you want to define QuAcc settings on-the-fly without writing them to a YAML file, you can instead modify the desired settings by defining individual environment variables with  `QUACC` as the prefix. For instance, to modify the `SCRATCH_DIR` setting to be `/tmp`, simply define `QUACC_SCRATCH_DIR=/tmp` as a new environment variable.
