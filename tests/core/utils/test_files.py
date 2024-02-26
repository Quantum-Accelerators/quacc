import gzip
import logging
import os
from pathlib import Path

import pytest

from quacc.utils.files import check_logfile, copy_decompress_files, make_unique_dir

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
@pytest.mark.parametrize("files_to_copy", ["src", ["src"], "sr*"])
def test_copy_decompress_files(tmp_path, files_to_copy):
    src = tmp_path / "src"
    src.mkdir()

    dst = tmp_path / "dst"
    dst.mkdir()

    Path(src / "file1").touch()
    Path(src / "dir1").mkdir()
    Path(src / "nested" / "nested").mkdir(parents=True)
    Path(src / "dir1" / "file2").touch()
    Path(src / "dir1" / "symlink1").symlink_to(src)

    copy_decompress_files(tmp_path, files_to_copy, dst)

    assert os.listdir(dst) == ["src"]
    assert sorted(os.listdir(dst / "src")) == ["dir1", "file1", "nested"]
    assert sorted(os.listdir(dst / "src" / "dir1")) == ["file2"]
    assert sorted(os.listdir(dst / "src" / "nested")) == ["nested"]


@pytest.mark.skipif(os.name == "nt", reason="Windows doesn't support symlinks")
@pytest.mark.parametrize("files_to_copy", ["file1", ["file1"], "file*"])
def test_copy_decompress_files_v2(tmp_path, files_to_copy):
    src = tmp_path / "src"
    src.mkdir()

    dst = tmp_path / "dst"
    dst.mkdir()

    Path(src / "file1").touch()
    Path(src / "dir1").mkdir()
    Path(src / "nested" / "nested").mkdir(parents=True)
    Path(src / "dir1" / "file2").touch()
    Path(src / "dir1" / "symlink1").symlink_to(src)

    copy_decompress_files(src, files_to_copy, dst)

    assert os.listdir(dst) == ["file1"]

    copy_decompress_files(src / "dir1", "file2", dst)

    assert sorted(os.listdir(dst)) == ["file1", "file2"]


@pytest.mark.skipif(os.name == "nt", reason="Windows doesn't support symlinks")
@pytest.mark.parametrize("files_to_copy", [Path("dir1", "file2"), "dir1/file2"])
def test_copy_decompress_files_v3(tmp_path, files_to_copy):
    src = tmp_path / "src"
    src.mkdir()

    dst = tmp_path / "dst"
    dst.mkdir()

    Path(src / "file1").touch()
    Path(src / "dir1").mkdir()
    Path(src / "nested" / "nested").mkdir(parents=True)
    Path(src / "dir1" / "file2").touch()
    Path(src / "dir1" / "symlink1").symlink_to(src)

    copy_decompress_files(src, files_to_copy, dst)

    assert os.listdir(dst) == ["dir1"]
    assert os.listdir(dst / "dir1") == ["file2"]


@pytest.mark.skipif(os.name == "nt", reason="Windows doesn't support symlinks")
def test_copy_decompress_files_v4(tmp_path):
    src = tmp_path / "src"
    src.mkdir()

    dst = tmp_path / "dst"
    dst.mkdir()

    Path(src / "file1").touch()
    Path(src / "dir1").mkdir()
    Path(src / "nested" / "nested").mkdir(parents=True)
    Path(src / "dir1" / "file2").touch()
    Path(src / "dir1" / "symlink1").symlink_to(src)

    copy_decompress_files(src / "dir1", "*", dst)

    assert os.listdir(dst) == ["file2"]


def test_copy_decompress_files_from_dir_warning(caplog):
    with caplog.at_level(logging.WARNING):
        copy_decompress_files("fake", "file", "test")
    assert "Cannot find file" in caplog.text


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

    copy_decompress_files(src, ["dir1/file2", "dir1/symlink1"], dst)

    assert (dst / "dir1" / "file2").exists()
    assert not (dst / "dir1" / "file1").exists()
    assert not (dst / "dir1" / "symlink1").exists()


def test_check_logfile(tmp_path):
    with open(tmp_path / "logs.out", "w") as f:
        f.write("trigger")
    assert check_logfile(str(tmp_path / "logs.out"), "trigger") is True
    assert check_logfile(str(tmp_path / "logs.out"), "test") is False

    with gzip.open(tmp_path / "logs2.out.gz", "wb") as gf:
        gf.write(b"trigger")
    assert check_logfile(str(tmp_path / "logs2.out"), "trigger") is True
    assert check_logfile(str(tmp_path / "logs2.out"), "test") is False
