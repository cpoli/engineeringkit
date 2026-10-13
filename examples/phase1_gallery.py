"""Small runnable Phase-1 gallery."""

import numpy as np
from scipy import signal

from engineeringkit.control_systems import PID, second_order_metrics
from engineeringkit.electrical_engineering import SeriesRLC, maximum_power_transfer
from engineeringkit.fluid_engineering import colebrook_friction_factor, reynolds_number
from engineeringkit.structural_engineering import PlanarTruss, SimplySupportedBeam
from engineeringkit.thermodynamics import OttoCycle, carnot_efficiency

beam = SimplySupportedBeam(5, 200e9, 8e-6).solve_point_load(10_000)
print("beam", beam.summary())
nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]])
truss = PlanarTruss(nodes, ((0, 1), (0, 2), (1, 2)), 200e9, 1e-3).solve([[0, 0], [0, 0], [0, -1000]], [0, 1, 3])
print("truss", truss.summary())
print("second order", second_order_metrics(0.5, 2))
PID(3, 1, 0.2).simulate(signal.TransferFunction([1], [1, 2, 1]))
SeriesRLC(10, 0.1, 100e-6).step_response()
print("max power", maximum_power_transfer(12, 6))
print("otto", OttoCycle(10).summary(), "carnot", carnot_efficiency(600, 300))
Re = reynolds_number(998, 2, 0.05, 1.002e-3)
print("pipe", Re, colebrook_friction_factor(Re, 1e-4))
