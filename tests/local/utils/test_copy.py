from pathlib import Path

from quacc.utils.files import copy_decompress


def test_copy_decompress(tmp_path):
    """Test copy_decompress function."""
    src = tmp_path / "src"
    src.mkdir()

    dst = tmp_path / "dst"
    dst.mkdir()

    Path(src / "file1").touch()
    Path(src / "dir1").mkdir()
    Path(f"{str(src)}{'/nested' * 10}").mkdir(parents=True)
    Path(src / "dir1" / "file2").touch()
    Path(src / "dir1" / "symlink1").symlink_to(src)

    copy_decompress(src, dst)

    assert (dst / "file1").exists()
    assert (dst / "dir1").exists()
    assert (dst / "dir1" / "file2").exists()
    assert Path(f"{str(dst)}{'/nested' * 10}").exists()
    assert not (dst / "dir1" / "symlink1").exists()
