# Running an automated testsuite with Python

## Prerequisites

Python 3.10 is recommended to ensure compatibility, other versions are probably 
compatible but not tested. Installation of the following packages is mandatory:

* [PyFoam](https://pypi.org/project/PyFoam/) 
* [pytest](https://pypi.org/project/pytest/)
* [caseFOAM](https://casefoam.readthedocs.io/en/latest/?badge=latest)
* [oftest](https://oftest.readthedocs.io/en/latest/)
* [pytest-xdist](https://pytest-xdist.readthedocs.io/en/latest/#pytest-xdist)

Typically installed with the command:

```bash
pip install -U <package-name>
```

## Usage

With a sourced working foam-extend 4.1, the test can be started by running the 
command below which executes the testing process from tests generation to results 
analysis. Check the "run" script for the individual steps.

```bash
./run
```

For cleaning up all test again run

```bash
./clean
```
