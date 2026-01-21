"""
Job argument classes for workflow engine interoperability.

This module provides a pattern for creating workflow-engine-specific
objects that can be used as arguments in jobs within a DAG.
The primary use case is copying files from source directories (which are
represented as futures in various workflow engines) to a target
directory.

Different workflow engines represent deferred/future values differently:
- Prefect/Parsl: Use `.result()` to resolve futures
- Dask: Uses `.compute()` to resolve delayed values
- Redun: Requires wrapping operations in `@task` decorators for lazy evaluation
- Jobflow: Uses UUID-based references that get resolved by the executor
- No engine (direct execution): Values are already resolved
"""

from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

from monty.json import MSONable

from quacc import get_settings
from quacc.utils.files import copy_decompress_files

if TYPE_CHECKING:
    from collections.abc import Mapping


class JobArgument:
    """
    Abstract factory base class for workflow-engine-specific job arguments.

    Subclasses must define a `_class_mapping` class variable that maps workflow
    engine names (str) to their corresponding implementation classes. The key
    `None` should map to the default implementation used when no workflow
    engine is configured.
    """

    def __new__(cls, *args, **kwargs):
        settings = get_settings()

        target_cls = cls._class_mapping.get(
            settings.WORKFLOW_ENGINE, cls._class_mapping[None]
        )

        obj = target_cls.__new__(target_cls)
        obj.__init__(*args, **kwargs)
        return obj


class JobflowCopy(MSONable):
    """
    Copy implementation for Jobflow engine.

    Parameters
    ----------
    futures_to_vals
        A mapping from Jobflow OutputReferences to the list of files to copy
        from the directory that the OutputReference will resolve to.
        Used when constructing the object before serialization.
    src_dir_to_files
        A mapping from resolved directory paths to file lists.
        Used when reconstructing the object after deserialization.
    """

    def __init__(self, futures_to_vals=None, src_dir_to_files=None):
        self.futures_to_vals = futures_to_vals
        self.src_dir_to_files = src_dir_to_files

    def as_dict(self):
        """
        Serialization method required for `MSONable` objects.

        The OutputReferences are stored in `uuids_to_refs`. Jobflow will
        resolve them before `from_dict` is called.
        """
        uuids_to_refs = {k.uuid: k for k in self.futures_to_vals}
        uuids_to_vals = {k.uuid: v for k, v in self.futures_to_vals.items()}
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "uuids_to_refs": uuids_to_refs,
            "uuids_to_vals": uuids_to_vals,
        }

    @classmethod
    def from_dict(cls, d):
        """
        Reconstruct from serialized dict after OutputReferences are resolved.
        """
        uuids_to_refs = d["uuids_to_refs"]
        uuids_to_vals = d["uuids_to_vals"]
        src_dir_to_files = {
            uuids_to_refs[uuid]: files for uuid, files in uuids_to_vals.items()
        }
        return cls(src_dir_to_files=src_dir_to_files)

    def do_copy(self, tgt_dir):
        """
        Copy files from all source directories to the target directory.

        Parameters
        ----------
        tgt_dir
            The destination directory where files will be copied.
        """
        for src_dir, files in self.src_dir_to_files.items():
            copy_decompress_files(src_dir, files, tgt_dir)


class DictCopy(dict):
    """
    Default copy implementation when no workflow engine is configured.

    This class inherits from dict so it can be used directly as a mapping of
    source_dir -> files_to_copy.
    """

    def do_copy(self, tgt_dir):
        """
        Copy files from all source directories to the target directory.

        Parameters
        ----------
        tgt_dir
            The destination directory where files will be copied.
        """
        for src_dir, files in self.items():
            copy_decompress_files(src_dir, files, tgt_dir)


class PrefectOrParslCopy:
    """
    Copy implementation for Prefect and Parsl workflow engines.

    Both Prefect and Parsl use futures that are resolved via the `.result()`
    method. This class stores a mapping from future objects to files, and
    resolves each future when `do_copy` is called.

    Parameters
    ----------
    src_dir_to_files
        A mapping from future objects (representing source directories) to
        lists of files to copy from those directories.
    """

    def __init__(self, src_dir_to_files):
        self.src_dir_to_files = src_dir_to_files

    def do_copy(self, tgt_dir):
        """
        Resolve futures and copy files to the target directory.

        Parameters
        ----------
        tgt_dir
            The destination directory where files will be copied.
        """
        for future, files in self.src_dir_to_files.items():
            copy_decompress_files(future.result(), files, tgt_dir)


class RedunCopy:
    """
    Copy implementation for the Redun workflow engine.

    Redun uses a lazy evaluation model where deferred values are only resolved
    when they cross a task boundary. This class wraps the copy operation in a
    redun `@task` to ensure that source directory values are properly resolved
    before the copy occurs.

    Parameters
    ----------
    src_dir_to_files
        A mapping from redun Value objects (representing source directories)
        to lists of files to copy from those directories.
    """

    def __init__(self, src_dir_to_files):
        self.src_dir_to_files = src_dir_to_files

    def do_copy(self, tgt_dir):
        """
        Create redun tasks to copy files to the target directory.

        Unlike other implementations, this returns a list of redun task results
        rather than performing the copy immediately. Redun's scheduler will
        execute these tasks when their input dependencies are resolved.

        Parameters
        ----------
        tgt_dir
            The destination directory where files will be copied.

        Returns
        -------
        list
            A list of redun task results (one per source directory).
        """
        from redun import task

        @task
        def _do_copy_task(src_dir, files, tgt_dir):
            # At this task boundary, redun evaluates src_dir from a deferred
            # value to its actual resolved value (the directory path string).
            copy_decompress_files(src_dir, files, tgt_dir)

        # Return a list of task invocations; redun will schedule these and
        # resolve the src_dir values when their upstream dependencies complete.
        return [
            _do_copy_task(src_dir, files, tgt_dir)
            for src_dir, files in self.src_dir_to_files.items()
        ]


class DaskCopy(dict):
    """
    Copy implementation for the Dask workflow engine.

    Dask uses delayed objects that are resolved via the `.compute()` method.
    This class inherits from dict to ensure proper serialization when passed
    between Dask workers.

    Parameters
    ----------
    src_dir_to_files
        A mapping from Dask delayed objects (representing source directories)
        to lists of files to copy from those directories.
    """

    def __init__(self, src_dir_to_files):
        # Store in the dict under a private key for serialization compatibility
        self["_src_dir_to_files"] = src_dir_to_files

    def do_copy(self, tgt_dir):
        """
        Compute delayed values and copy files to the target directory.

        Parameters
        ----------
        tgt_dir
            The destination directory where files will be copied.
        """
        for future, files in self["_src_dir_to_files"].items():
            copy_decompress_files(future.compute(), files, tgt_dir)


class Copy(JobArgument):
    """
    Unified API for copying files from job output directories in a workflow DAG.

    This is the main class that users should instantiate. It acts as a factory
    that automatically returns the appropriate engine-specific implementation
    based on the current WORKFLOW_ENGINE setting.
    """

    # Maps workflow engine names to their implementation classes.
    # The None key provides the default for when no engine is configured.
    _class_mapping: ClassVar[Mapping[str | None, type[JobArgument]]] = {
        "jobflow": JobflowCopy,
        "prefect": PrefectOrParslCopy,
        "parsl": PrefectOrParslCopy,
        "redun": RedunCopy,
        "dask": DaskCopy,
        None: DictCopy,
    }
