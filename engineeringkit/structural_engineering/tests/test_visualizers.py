"""Smoke tests for engineeringkit.structural_engineering.visualizers and the result plots."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from engineeringkit import structural_engineering as s
from engineeringkit.structural_engineering.visualizers import plot_beam_diagrams, plot_frame, plot_mohr_circle, plot_truss


def test_beam_diagrams_need_three_axes():
    result = s.SimplySupportedBeam(4.0, 200e9, 1e-5).solve_uniform_load(1e3)
    _, axes = plt.subplots(2, 1)
    with pytest.raises(ValueError):
        plot_beam_diagrams(result, axes)
    assert len(plot_beam_diagrams(result)) == 3
    plt.close("all")


def test_truss_colours_members_by_force():
    nodes = [[0, 0], [1, 0], [0.5, 1]]
    ax = plot_truss(nodes, ((0, 1), (0, 2), (1, 2)), forces=[1.0, -1.0, 0.0])
    assert [line.get_color() for line in ax.lines] == ["C3", "C0", "k"]
    with pytest.raises(ValueError):
        plot_truss([0, 1], ((0, 1),))
    plt.close("all")


def test_frame_and_moment_distribution_plots():
    frame = s.PlanarFrame([[0, 0], [0, 3], [4, 3]], ((0, 1), (1, 2)), 200e9, 1e-2, 1e-4)
    result = frame.solve([0, 0, 0, 1e3, 0, 0, 0, -1e3, 0], [0, 1, 2])
    fig, ax = result.plot(scale=50)
    assert len(ax.lines) == 4
    assert plot_frame(result, ax=ax) is ax
    md = s.moment_distribution([(0, 1, 5.0, 1e7, 1e3), (1, 2, 5.0, 1e7, 1e3)], fixed_joints=(0, 2))
    fig, ax = md.plot()
    assert ax.get_yscale() == "log"
    plt.close("all")


def test_mohr_circle_marks_faces_and_principal_stresses():
    ax = plot_mohr_circle(80e6, -20e6, 30e6)
    assert ax.yaxis_inverted()
    assert len(ax.get_legend().get_texts()) == 2
    plt.close("all")
