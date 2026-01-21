from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

from quacc import get_settings
from quacc.utils.files import copy_decompress_files

if TYPE_CHECKING:
    from collections.abc import Mapping


class JobArgument:
    def __new__(cls, *args, **kwargs):
        settings = get_settings()

        target_cls = cls._class_mapping.get(
            settings.WORKFLOW_ENGINE, cls._class_mapping[None]
        )

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

        return [
            _do_copy_task(src_dir, files, tgt_dir)
            for src_dir, files in self.src_dir_to_files.items()
        ]


class DaskCopy(dict):
    def __init__(self, src_dir_to_files):
        self["_src_dir_to_files"] = src_dir_to_files

    def do_copy(self, tgt_dir):
        for future, files in self["_src_dir_to_files"].items():
            copy_decompress_files(future.compute(), files, tgt_dir)


class Copy(JobArgument):
    _class_mapping: ClassVar[Mapping[str | None, type[JobArgument]]] = {
        "jobflow": JobflowCopy,
        "prefect": PrefectOrParslCopy,
        "parsl": PrefectOrParslCopy,
        "redun": RedunCopy,
        "dask": DaskCopy,
        None: DictCopy,
    }
