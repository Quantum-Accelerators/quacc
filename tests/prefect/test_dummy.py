from prefect.testing.utilities import prefect_test_harness
from prefect import flow


@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    with prefect_test_harness():
        yield

def test_dummy():

    @flow
    def my_favorite_function():
        return 42

    assert my_favorite_function() == 42
