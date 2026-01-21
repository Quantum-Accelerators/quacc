from __future__ import annotations

from monty.json import MSONable

from quacc import get_settings
from quacc.utils.files import copy_decompress_files


class JobArgument:
    def __new__(cls, *args, **kwargs):
        settings = get_settings()

        match settings.WORKFLOW_ENGINE:
            case "jobflow":
                target_cls = cls._jobflow_cls
            case "prefect":
                target_cls = cls._prefect_cls
            case "parsl":
                target_cls = cls._parsl_cls
            case "redun":
                target_cls = cls._redun_cls
            case "dask":
                target_cls = cls._dask_cls
            case _:
                target_cls = cls._default_cls

        obj = target_cls.__new__(target_cls)
        obj.__init__(*args, **kwargs)
        return obj


class JobflowCopy:
    def __init__(self, futures_to_vals):
        self.futures_to_vals = futures_to_vals

    def as_dict(self):
        uuids_to_keys = {k.uuid: k for k in self.futures_to_vals}
        uuids_to_vals = {k.uuid: v for k, v in self.futures_to_vals.items()}
        return {"uuids_to_keys": uuids_to_keys, "uuids_to_vals": uuids_to_vals}


class DictCopy(dict):
    def do_copy(self, tgt_dir):
        for src_dir, files in self.items():
            copy_decompress_files(src_dir, files, tgt_dir)


class PrefectCopy:
    def __init__(self, src_dir_to_files):
        self.src_dir_to_files = src_dir_to_files

    def do_copy(self, tgt_dir):
        for src_dir, files in self.items():
            copy_decompress_files(src_dir, files, tgt_dir)


class PrefectOrParslCopy:
    def __init__(self, src_dir_to_files):
        self.src_dir_to_files = src_dir_to_files

    def do_copy(self, tgt_dir):
        for future, files in self.src_dir_to_files.items():
            copy_decompress_files(future.result(), files, tgt_dir)


class RedunCopy:
    def __init__(self, src_dir_to_files):
        self.src_dir_to_files = src_dir_to_files

    def do_copy(self, tgt_dir):
        from redun import task

        @task
        def _do_copy_task(src_dir, files, tgt_dir):
            # At THIS boundary, redun evaluates src_dir
            copy_decompress_files(src_dir, files, tgt_dir)

        tasks = [_do_copy_task(src_dir, files, tgt_dir) for src_dir, files in self.src_dir_to_files.items()]
        return tasks


class DaskCopy(dict):
    def __init__(self, src_dir_to_files):
        self["_src_dir_to_files"] = src_dir_to_files

    def do_copy(self, tgt_dir):
        for future, files in self['_src_dir_to_files'].items():
            copy_decompress_files(future.compute(), files, tgt_dir)


class Copy(JobArgument):
    _jobflow_cls = JobflowCopy
    _prefect_cls = PrefectOrParslCopy
    _parsl_cls = PrefectOrParslCopy
    _redun_cls = RedunCopy
    _dask_cls = DaskCopy
    _default_cls = DictCopy
