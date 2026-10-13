"""Compact executable gallery touching every EngineeringKit domain."""

from engineeringkit.aerodynamics import thin_airfoil_cl
from engineeringkit.communications import shannon_capacity
from engineeringkit.heat_transfer import fourier_heat_rate
from engineeringkit.materials import critical_crack_length
from engineeringkit.mechanical_engineering import natural_frequency_hz
from engineeringkit.reliability import series_reliability
from engineeringkit.signal_processing import nyquist_frequency

print("heat W", fourier_heat_rate(15, 1, 20, 0.1))
print("critical crack m", critical_crack_length(50e6, 100e6))
print("Nyquist Hz", nyquist_frequency(1000))
print("natural Hz", natural_frequency_hz(1000, 10))
print("CL", thin_airfoil_cl(0.05))
print("capacity bit/s", shannon_capacity(1e6, 10))
print("series reliability", series_reliability(0.99, 0.98, 0.97))
