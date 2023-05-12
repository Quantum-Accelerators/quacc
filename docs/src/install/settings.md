# Settings

## Introduction

If you've completed the prior instructions, you're done and can run QuAcc! As you use QuAcc, there may be a few default settings you wish to change. For instance, you can see that by default, `GZIP_FILES` is set to `True`, which means that QuAcc will automatically compress the output files from the calculation. Similarly, you can see that `SCRATCH_DIR` is set to the `$SCRATCH` environment variable if it exists.

## Modifying the Default Settings

By defualt, the global QuAcc settings can be found in the `settings.py` file [here](https://github.com/arosen93/quacc/blob/main/quacc/settings.py). If you wish to modify any of the global QuAcc settings, there are two ways to do so.

### Method 1

The easiest way is to create a custom YAML file with custom settings. For instance, it might look something like the following to modify the aforementioned settings.

```yaml
SCRATCH_DIR: /tmp
GZIP_FILES: false
```

By default, QuAcc looks for this YAML file at `~/quacc.yaml`. If you wish to use astore the `quacc.yaml` file somewhere else, you can define the environment variable `QUACC_CONFIG_FILE` and point it to the `quacc.yaml` path of your choosing.

### Method 2

If you want to define QuAcc settings on-the-fly without writing them to a YAML file, you can instead modify the desired settings by defining individual environment variables with  `QUACC` as the prefix. For instance, to modify the `SCRATCH_DIR` setting to be `/tmp`, simply define `QUACC_SCRATCH_DIR=/tmp` as a new environment variable.
