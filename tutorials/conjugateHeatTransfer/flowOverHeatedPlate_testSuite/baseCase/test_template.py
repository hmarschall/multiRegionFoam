import os
import pytest
import oftest
from oftest import run_reset_case

def test_completed(run_reset_case):
    log = oftest.path_log()
    assert oftest.case_status(log) == "completed"
