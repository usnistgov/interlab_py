"""Tests for `interlab` package."""

from __future__ import annotations


def test_version() -> None:
    from interlab import __version__

    assert __version__ != "999"
