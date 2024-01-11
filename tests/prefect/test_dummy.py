from prefect import flow
from prefect.testing.utilities import prefect_test_harness


@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    with prefect_test_harness():
        yield


def test_dummy():
    @flow
    def my_favorite_function():
        return 42

    assert my_favorite_function() == 42
