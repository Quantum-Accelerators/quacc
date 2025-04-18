# Results Handling

It can be useful to understand how quacc handles files so that you can adjust the behavior for your needs. Let's consider a recipe that writes an input file `INPUT` and an output file `OUTPUT`.

!!! Tip "Having an Organized Filesystem"

    If you want to achieve a bit more structure in how your files are stored, you can dynamically set the `RESULTS_DIR` for each project of interest. Importantly, this should be done by modifying the `QUACC_RESULTS_DIR` environment variable (e.g. in your submission script) rather than modifying the YAML file directly, such that the change reflects only your given workflow. Similarly, to prevent the automatic creation of a unique directory for each job, you can set the `QUACC_CREATE_UNIQUE_DIR` environment variable to `False`.

## Scenario 1: The Default Case

Here, let us begin by assuming we are using the default quacc settings. By default, `RESULTS_DIR` is set to `Path.cwd()` (i.e. the current working directory), but this can be modified by the user.

### Job Runtime

At job runtime, the file structure looks like:

```text
RESULTS_DIR
├── tmp-quacc-2023-12-08-67890
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

!!! Tip "Disabling Automated File Compression"

    If you don't want the results to be gzipped, you can set the `GZIP_FILES` setting to `False`.

### Job Failure

If the job fails or does not complete, then the file structure looks like:

```text
RESULTS_DIR
├── failed-quacc-2023-12-08-67890
│   ├── INPUT
    └── OUTPUT
```

## Scenario 2: Specifying a `SCRATCH_DIR`

By default, the `SCRATCH_DIR` setting is set to `None`. This means that all the calculations will be run somewhere within the `RESULTS_DIR`. Sometimes, the filesystem where you want to permanently store the results and where you want to run the calculations is different. In this case, you can specify a `SCRATCH_DIR` setting to distinguish between the two.

Here, let's assume the user has specified the `SCRATCH_DIR` setting to be a custom path.

### Job Runtime

At job runtime, the file structure looks like:

```text
RESULTS_DIR
├── symlink-tmp-quacc-2023-12-08-67890
│
```

```text
SCRATCH_DIR
├── tmp-quacc-2023-12-08-67890
│   ├── INPUT
    └── OUTPUT
```

Here, the `symlink-tmp-quacc-2023-12-08-67890` is a temporary symbolic link that points to `SCRATCH_DIR/tmp-quacc-2023-12-08-67890` so you can easily monitor the progress of the calculation. On Windows, no symbolink link is created, but the `tmp-quacc-2023-12-08-67890` directory is still created in `SCRATCH_DIR`.

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

If the job fails or does not complete, then the file structure looks like:

```text
RESULTS_DIR
├── symlink-failed-quacc-2023-12-08-67890
│
```

```text
SCRATCH_DIR
├── failed-quacc-2023-12-08-67890
│   ├── INPUT
    └── OUTPUT
```
