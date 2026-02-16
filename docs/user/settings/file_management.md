# File Management

It can be useful to understand how quacc handles files so that you can adjust the behavior for your needs. Let's consider a recipe that writes an input file `INPUT` and an output file `OUTPUT`.

!!! Tip "Having an Organized Filesystem"

    If you want to achieve a bit more structure in how your files are stored, you can dynamically set the `RESULTS_DIR` quacc setting as you wish. If you wish to prevent the automatic creation of a unique directory for each job, you can set the `CREATE_UNIQUE_DIR` quacc setting to `False`.

    For a hierarchical folder structure that reflects the relationship between flows, subflows, and jobs, see [Auto-generation of folder structure](#scenario-2-auto-generation-of-folder-structure).

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

## Scenario 2: Auto-generation of folder structure

Quacc has a new feature for "auto-generation" of a sensible output folder structure during job execution. To enable this feature, set `AUTODISCOVER_DIR` to `true`.

Here, let us assume that `RESULTS_DIR` is set to `Path.cwd()` (i.e. the current working directory), though this can be modified just like any other setting.

### Job Success

Once the job successfully completes, the file structure looks like:

```text
RESULTS_DIR
├── quacc-2026-02-16-67890
    └── relax_job-2026-02-16-71530
        ├── INPUT.gz
        └── OUTPUT.gz
```

Here, `relax_job` is the name of the job and is suffixed with a timestamp. A recipe may have multiple jobs with the same name, and the timestamp helps us distinguish one run of a job with another.

For typical quacc recipes that have a flow, subflow(s), and job(s), the output folders will reflect the relationship between them. For example, the `bulk_to_slabs_flow` recipe results in something like:

```text
    RESULTS_DIR
    └── quacc-2026-02-16-25313
        └── bulk_to_slabs_flow-2026-02-16-35312
            └── bulk_to_slabs_subflow-2026-02-16-38381
                ├── relax_job-2026-02-16-45812
                │         ├── opt.json.gz
                │         ├── opt.log.gz
                │         └── opt.traj.gz
                ├── relax_job-2026-02-16-42525
                │         ├── opt.json.gz
                │         ├── opt.log.gz
                │         └── opt.traj.gz
                ...
```

At each level, the folder name has a prefix that identifies the name of the `@flow`, `@subflow`, or `@job` being run, and the suffix indicates the timestamp. In general, the ordering of folders at a level corresponds to the order in which they were run.

## Specifying a `SCRATCH_DIR`

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

If `AUTODISCOVER_DIR` is `false`:

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

If `AUTODISCOVER_DIR` is `true`:

```text
RESULTS_DIR
├── quacc-2026-02-16-67890
    └── relax_job-2026-02-16-71530
        ├── INPUT.gz
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
