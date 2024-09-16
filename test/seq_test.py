import pytest
from SOM_Seq_Sim.Seq_Sim.utils import *


def test_optimal_layout():
    assert optimal_subplot_layout(9) == (3, 3)
    assert optimal_subplot_layout(10) == (2, 5)
    assert optimal_subplot_layout(11) == (6, 2)
    assert optimal_subplot_layout(12) == (2, 6)


def test_remainder():
    assert remainder(10, 3) == 1
    assert remainder(12, 5) == 2
    assert remainder(7, 2) == 1
    assert remainder(6, 3) == 0
