# File Management

It can be useful to understand how quacc handles files so that you can adjust the behavior for your needs. Let's consider a recipe that writes an input file INPUT and an output file OUTPUT.

!!! Hint "Having an Organized Filesystem"

    If you want to achieve a bit more structure in how your files are stored, you can dynamically set the `RESULTS_DIR` for each project of interest. Importantly, this should be done by modifying the `QUACC_RESULTS_DIR` environment variable (e.g. in your submission script) rather than modifying the YAML file directly, such that the change reflects only your given workflow.

## Scenario 1: The Default Case

Here, let us begin by assuming we are using the default quacc settings. By default, `RESULTS_DIR` is set to `Path.cwd()` (i.e. the current working directory), but this can be modified by the user.

### Job Runtime

At job runtime, the file structure looks like:

```text
RESULTS_DIR
├── quacc-tmp-12345
│   ├── INPUT
    └── OUTPUT
```

### Job Success

Once the job successfully completes, the file structure looks like:

```text
RESULTS_DIR
├── INPUT.gz
├── OUTPUT.gz
│
```

!!! Tip

    If you don't want the results to be gzipped, you can set the `GZIP_FILES` setting to `False`.

### Job Failure

If the job fails or does not complete, then the `quacc-tmp-12345` directory will remain in `RESULTS_DIR` so you can inspect the files.

## Scenario 2: Specifying a `SCRATCH_DIR`

By default, the `SCRATCH_DIR` setting is set to `None`. This means that all the calculations will be run somewhere within the `RESULTS_DIR`. Sometimes, the filesystem where you want to permanently store the results and where you want to run the calculations are different. In this case, you can specify a `SCRATCH_DIR` setting to distinguish between the two.

Here, let's assume the user has specified the `SCRATCH_DIR` setting to be a custom `Path`.

### Job Runtime

At job runtime, the file structure looks like:

```text
RESULTS_DIR
├── quacc-tmp-12345-symlink
│
```

```text
SCRATCH_DIR
├── quacc-tmp-12345
│   ├── INPUT
    └── OUTPUT
```

Here, the `quacc-tmp-12345-symlink` is a temporary symbolic link that points to `SCRATCH_DIR/quacc-tmp-12345` so you can easily monitor the progress of the calculation.

### Job Success

Once the job successfully completes, the file structure looks like:

```text
RESULTS_DIR
├── INPUT.gz
├── OUTPUT.gz
│
```

```text
SCRATCH_DIR
│
```

### Job Failure

If the job fails or does not complete, then the `quacc-tmp-12345` directory will remain in `SCRATCH_DIR` so you can inspect the files. The symbolic link in `RESULTS_DIR` will also remain.

## Scenario 3: Setting `CREATE_UNIQUE_DIR` to `True`

If you run multiple jobs, the resulting files will start to overwrite each other in `RESULTS_DIR` by default. Some workflow engines have settings to circumvent this, but regardless, you can set the `CREATE_UNIQUE_DIR` setting to `True` in quacc to ensure no file clashes occur.

Let's revisit Scenario 2 and assume the user has also set `CREATE_UNIQUE_DIR` to `True`.

### Job Runtime

At job runtime, the file structure looks like:

```text
RESULTS_DIR
├── quacc-2023-12-08-67890
│   └── quacc-tmp-12345-symlink
```

```text
SCRATCH_DIR
├── quacc-tmp-12345
│   ├── INPUT
    └── OUTPUT
```

### Job Success

Once the job successfully completes, the file structure looks like:

```text
RESULTS_DIR
├── quacc-2023-12-08-67890
│   ├── INPUT.gz
    └── OUTPUT.gz
```

```text
SCRATCH_DIR
│
```

### Job Failure

If the job fails or does not complete, then the `quacc-tmp-12345` directory will remain in `SCRATCH_DIR` so you can inspect the files. The symbolic link in `RESULTS_DIR/quacc-2023-12-08-67890` will also remain.
