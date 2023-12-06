from ase.build import bulk

from quacc import SETTINGS
from quacc.recipes.espresso.core import ph_job, static_job

DEFAULT_SETTINGS = SETTINGS.model_copy()


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk('Si')

    input_data = {
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.005,
        'mixing_mode': 'plain',
        'mixing_beta': 0.6,
    }

    pseudopotentials = {'Si': 'Si.upf'}
    
    results = static_job(atoms,
                        input_data=input_data,
                        pseudopotentials=pseudopotentials,
                        k_spacing=0.5)
    
    assert results
    
def test_ph_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk('Li')

    input_data = {
        'calculation': 'vc-relax',
        'occupations': 'smearing',
        'smearing': 'cold',
        'degauss': 0.02,
        'mixing_mode': 'TF',
        'mixing_beta': 0.7,
    }

    ph_loose = {
        'tr2_ph': 1e-10,
    }

    pseudopotentials = {'Li': 'Li.upf'}

    pw_results = static_job(atoms,
                        input_data=input_data,
                        pseudopotentials=pseudopotentials,
                        k_spacing=0.25)
    
    ph_results = ph_job(input_data=ph_loose,
                    copy_files=pw_results['dir_name'])
    
    assert ph_results