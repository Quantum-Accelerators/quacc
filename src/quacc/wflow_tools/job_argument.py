from __future__ import annotations

from monty.json import MSONable

from quacc import get_settings
from quacc.utils.files import copy_decompress_files


class JobArgument:
    def __new__(cls, **kwargs):
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
        obj.__init__(**kwargs)
        return obj


class JobflowCopy(MSONable):
    def __init__(self, src_dir, files: str = "*"):
        self.src_dir = src_dir
        self.files = files

    def do_copy(self, tgt_dir):
        copy_decompress_files(self.src_dir, self.files, tgt_dir)


class DictCopy(dict):
    def do_copy(self, tgt_dir):
        copy_decompress_files(self["src_dir"], self.get("files", "*"), tgt_dir)


class ParslCopy:
    def __init__(self, src_dir, files: str = "*"):
        self.src_dir = src_dir
        self.files = files

    def do_copy(self, tgt_dir):
        copy_decompress_files(self.src_dir.result(), self.files, tgt_dir)


class RedunCopy:
    def __init__(self, src_dir, files: str = "*"):
        self._src_dir = src_dir
        self._files = files

    def do_copy(self, tgt_dir):
        from redun import task

        @task
        def _do_copy_task(src_dir, files, tgt_dir):
            # At THIS boundary, redun evaluates src_dir
            copy_decompress_files(src_dir, files, tgt_dir)

        return _do_copy_task(self._src_dir, self._files, tgt_dir)


class DaskCopy(dict):
    def __init__(self, src_dir, files: str = "*"):
        self["_src_dir"] = src_dir
        self["_files"] = files

    def do_copy(self, tgt_dir):
        copy_decompress_files(self["_src_dir"].compute(), self["_files"], tgt_dir)


class Copy(JobArgument):
    _jobflow_cls = JobflowCopy
    _prefect_cls = DictCopy
    _parsl_cls = ParslCopy
    _redun_cls = RedunCopy
    _dask_cls = DaskCopy
    _default_cls = DictCopy
