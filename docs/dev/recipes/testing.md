# Testing Recipes

All recipes should be tested by adding a corresponding suite of tests to `quacc/tests/local/recipes`.

The "local" tests do not actually use any workflow orchestration, but this is often a good place to add thorough unit tests for various features of the code.

Additionally, if you use a workflow engine, it is strongly advised to add a corresponding unit test, e.g. in `quacc/tests/parsl` if you use Parsl. This test can be a simple functional test rather than a detailed suite of unit tests.