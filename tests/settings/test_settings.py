from pathlib import Path

from quacc import SETTINGS

DEFAULT_SETTINGS = SETTINGS.copy()


FILE_DIR = Path(__file__).resolve().parent


def setup_module():
    SETTINGS.WORKFLOW_ENGINE = "local"


def teardown_module():
    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE


def setup_function():
    SETTINGS.PRIMARY_STORE = None
    SETTINGS.GZIP_FILES = True
    SETTINGS.CREATE_UNIQUE_WORKDIR = False


def teardown_function():
    SETTINGS.PRIMARY_STORE = DEFAULT_SETTINGS.PRIMARY_STORE
    SETTINGS.GZIP_FILES = DEFAULT_SETTINGS.GZIP_FILES
    SETTINGS.CREATE_UNIQUE_WORKDIR = DEFAULT_SETTINGS.CREATE_UNIQUE_WORKDIR


def test_file(monkeypatch, tmpdir):
    import os

    from quacc.settings import QuaccSettings

    tmpdir.chdir()

    assert QuaccSettings().GZIP_FILES is True

    with open("quacc_test.yaml", "w") as f:
        f.write("GZIP_FILES: false")
    monkeypatch.setenv(
        "QUACC_CONFIG_FILE", os.path.join(os.getcwd(), "quacc_test.yaml")
    )

    assert QuaccSettings().GZIP_FILES is False
    os.remove("quacc_test.yaml")


def test_store(tmpdir):
    from ase.build import bulk
    from maggma.stores import MemoryStore

    from quacc.recipes.emt.core import static_job

    tmpdir.chdir()
    SETTINGS.PRIMARY_STORE = MemoryStore()
    atoms = bulk("Cu")
    static_job(atoms)


def test_results_dir(tmpdir):
    import os

    from ase.build import bulk

    from quacc.recipes.emt.core import relax_job

    tmpdir.chdir()

    atoms = bulk("Cu")
    relax_job(atoms)
    assert "opt.traj.gz" in os.listdir(os.getcwd())
    os.remove("opt.traj.gz")
    SETTINGS.GZIP_FILES = False
    relax_job(atoms)
    assert "opt.traj" in os.listdir(os.getcwd())
    os.remove("opt.traj")


def test_env_var(monkeypatch):
    from quacc.settings import QuaccSettings

    p = FILE_DIR / "my/scratch/dir"
    monkeypatch.setenv("QUACC_SCRATCH_DIR", p)
    assert p.expanduser().resolve() == QuaccSettings().SCRATCH_DIR


def test_yaml(tmpdir, monkeypatch):
    from quacc.settings import QuaccSettings

    tmpdir.chdir()

    p = FILE_DIR / "my/new/scratch/dir"
    with open("quacc_test.yaml", "w") as f:
        f.write(f"SCRATCH_DIR: {p}")
    monkeypatch.setenv("QUACC_CONFIG_FILE", "quacc_test.yaml")
    assert p.expanduser().resolve() == QuaccSettings().SCRATCH_DIR
