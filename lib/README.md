# renewablepotentialslib

Library code of the solar-and-wind-potentials workflow.

The library code contains general-purpose functions and routines that we expect to change at a slow pace -- in contrast to scripts. If you alter library code be aware that this will not trigger reruns of workflow rules. Think of `renewablepotentialslib` as any other dependency of this workflow (like NumPy): changing the version of any dependency will not rerun worfklow rules. When you change library code, you will have to rerun rules manually where needed.

## Developer Guide

### Installation

Best install `renewablepotentialslib` in editable mode:

    $ conda env create -f requirements-test.yaml
    $ conda activate renewablepotentialslib-test

### Run the test suite

Run the test suite with py.test:

    $ py.test
