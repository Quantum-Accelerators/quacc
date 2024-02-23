import gzip
import logging
import os
from pathlib import Path

import pytest

from quacc.utils.files import (
    check_logfile,
    copy_decompress_files_from_dir,
    copy_decompress_tree,
    make_unique_dir,
)

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = True


def test_make_unique_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    jobdir = make_unique_dir()
    assert os.path.exists(jobdir)

    jobdir = make_unique_dir(base_path="tmp_dir")
    assert os.path.exists("tmp_dir")
    assert "tmp_dir" in str(jobdir)
    assert os.path.exists(jobdir)


@pytest.mark.skipif(os.name == "nt", reason="Windows doesn't support symlinks")
def test_copy_decompress_files_from_dir(tmp_path):
    src = tmp_path / "src"
    src.mkdir()

    dst = tmp_path / "dst"
    dst.mkdir()

    Path(src / "file1").touch()
    Path(src / "dir1").mkdir()
    Path(f"{src}{'/nested' * 10}").mkdir(parents=True)
    Path(src / "dir1" / "file2").touch()
    Path(src / "dir1" / "symlink1").symlink_to(src)

    copy_decompress_files_from_dir(src, dst)

    assert (dst / "file1").exists()
    assert (dst / "dir1").exists()
    assert (dst / "dir1" / "file2").exists()
    assert Path(f"{dst}{'/nested' * 10}").exists()
    assert not (dst / "dir1" / "symlink1").exists()


@pytest.mark.skipif(os.name == "nt", reason="Windows doesn't support symlinks")
def test_copy_decompress_tree(tmp_path):
    src = tmp_path / "src"
    src.mkdir()

    dst = tmp_path / "dst"
    dst.mkdir()

    Path(src / "file1").touch()
    Path(src / "dir1").mkdir()
    Path(src / "dir1" / "dir2").mkdir()
    Path(src / "dir1" / "file2").touch()
    Path(src / "dir1" / "file1").touch()
    Path(src / "dir1" / "symlink1").symlink_to(src)

    to_copy = {Path(src): ["dir1/file2", "dir1/symlink1"]}

    copy_decompress_tree(to_copy, dst)

    assert (dst / "dir1" / "file2").exists()
    assert not (dst / "dir1" / "file1").exists()
    assert not (dst / "dir1" / "symlink1").exists()

    to_copy = {Path(src, "dir1", "dir2"): ["../file2", "../symlink1"]}

    copy_decompress_tree(to_copy, dst / "dir1")

    assert Path(dst, "file2").exists()
    assert not Path(dst, "file1").exists()


@pytest.mark.skipif(os.name == "nt", reason="Windows doesn't support symlinks")
def test_copy_decompress_files_from_dir_v2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    src = Path(os.getcwd())

    dst = tmp_path / "dst"
    dst.mkdir()

    Path(src / "file1").touch()
    Path(src / "dir1").mkdir()
    Path(f"{src}{'/nested' * 10}").mkdir(parents=True)
    Path(src / "dir1" / "file2").touch()
    Path(src / "dir1" / "symlink1").symlink_to(src)

    copy_decompress_files_from_dir(src, dst)

    assert (dst / "file1").exists()
    assert (dst / "dir1").exists()
    assert (dst / "dir1" / "file2").exists()
    assert Path(f"{dst}{'/nested' * 10}").exists()
    assert not (dst / "dir1" / "symlink1").exists()


def test_copy_decompress_files_from_dir_v3(caplog):
    with caplog.at_level(logging.WARNING):
        copy_decompress_files_from_dir("fake", "test")


def test_check_logfile(tmp_path):
    with open(tmp_path / "logs.out", "w") as f:
        f.write("trigger")
    assert check_logfile(str(tmp_path / "logs.out"), "trigger") is True
    assert check_logfile(str(tmp_path / "logs.out"), "test") is False

    with gzip.open(tmp_path / "logs2.out.gz", "wb") as gf:
        gf.write(b"trigger")
    assert check_logfile(str(tmp_path / "logs2.out"), "trigger") is True
    assert check_logfile(str(tmp_path / "logs2.out"), "test") is False
