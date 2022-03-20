# Settings

## Introduction

If you've completed the prior instructions, you're done and can run QuAcc! As you use QuAcc, there may be a few default settings you wish to change. We discuss how to do so here.

By defualt, the global QuAcc settings can be found in the `settings.py` file [here](https://github.com/arosen93/quacc/blob/main/quacc/settings.py). For instance, you can see that by default, `GZIP_FILES` is set to `True`, which means that QuAcc will automatically compress the output files from the calculation. Similarly, you can see that `SCRATCH_DIR` is set to the `$SCRATCH` environment variable if it exists.

## Modifying the Default Settings

If you wish to modify any of the QuAcc settings, there are several ways to do so.

The easiest way is to create a file `~/quacc.yaml`. For instance, it might looke something like the following to modify the aforementioned settings.

```yaml
SCRATCH_DIR: /tmp
GZIP_FILES: false
```

If you wish to store the `quacc.yaml` file somewhere else, you can define the environment variable `QUACC_CONFIG_FILE` and point it to the `quacc.yaml` path of your choosing.

Finally, you can also modify the settings by defining individual environment variables with the `QUACC` prefix (e.g. `QUACC_SCRATCH_DIR=/tmp`).
