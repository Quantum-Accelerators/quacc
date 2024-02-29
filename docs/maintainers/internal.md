# Overview

## Approving Jenkins Builds

We have a [Jenkins build pipeline](https://jenkins.princeton.edu/job/rosen_group/) that is automatically launched on pull requests, running a subset of the tests on the [Adroit cluster](https://researchcomputing.princeton.edu/systems/adroit) at Princeton with the `--noconftest` pytest flag. This makes it possible to test recipes for licensed codes in a true production environment.

When a new pull request is opened by someone who is on the whitelist (by default, everyone in the [@Quantum-Accelerators](https://github.com/Quantum-Accelerators) organization), the build process will automatically start.

For a pull request opened by a contributor not on the whitelist, [@buildbot-princeton](https://github.com/buildbot-princeton) will ask "Can one of the admins verify this patch?". No build will be launched until one of the following actions are taken by a member on the whitelist:

- "@buildbot-princeton test this please" to run the test once. Use this if you don't personally know the contributor. Use "@buildbot-princeton retest this please" to launch a follow-up test, if needed.
- "@buildbot-princeton ok to test" to run tests for all future commits in the pull request. Use this for trusted contributors only.
- "@buildbot-princeton add to whitelist" to add a contributor to the whitelist permanently. Use this only after consulting with a maintainer of `quacc`.

## Manual Upload to PyPI

To upload to PyPI manually in the case of a GitHub Actions failure:

```bash
pip install -e .
python -m build
twine check dist/*
twine upload dist/*
```
