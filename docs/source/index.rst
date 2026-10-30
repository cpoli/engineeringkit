engineeringkit: Computational Engineering in Python
===================================================

.. include:: /_generated/vars.rst

**engineeringkit** makes the models behind |num_subpackages| engineering domains
computable, inspectable and visualizable. Choose a domain to find its foundational
ideas, runnable examples and full API in one place. NumPy and SciPy provide mature
numerical solvers; EngineeringKit exposes engineering equations, SI conventions
and intermediate results for learning and teaching.

.. image:: _static/images/readme_hero.png
   :alt: Panel-method streamlines over a NACA 2412, two chaotic double pendulums released 0.001 rad apart, and the Van der Pol limit cycle computed with EngineeringKit
   :width: 100%

Each domain at a glance
-----------------------

- :mod:`engineeringkit.aerodynamics` -- Airfoil and wing baselines, drag polars and compressible flow.

  .. image:: _static/images/readme_aerodynamics.png
     :alt: Panel-method pressure distribution on a NACA 2412, lifting-line spanwise loading, and the oblique-shock theta-beta-Mach chart
     :width: 100%

- :mod:`engineeringkit.communications` -- Modulation and pulse shaping, noise, capacity, detection, PCM, source and channel coding, fading and link budgets.

  .. image:: _static/images/readme_communications.png
     :alt: A noisy 16-QAM constellation, a raised-cosine eye diagram, and bit error rate curves
     :width: 100%

- :mod:`engineeringkit.control_systems` -- Classical stability, compensator design, PID with anti-windup, nonlinear loops, state feedback, LQR and estimation.

  .. image:: _static/images/readme_control_systems.png
     :alt: PID step responses, a root locus, and a Nyquist plot
     :width: 100%

- :mod:`engineeringkit.electrical_engineering` -- DC and AC circuits, transients, electronics, machines and power systems.

  .. image:: _static/images/readme_electrical_engineering.png
     :alt: Series RLC step responses, a Smith chart, and induction motor torque-speed curves
     :width: 100%

- :mod:`engineeringkit.fluid_engineering` -- Hydrostatics, flow measurement, momentum, pipe friction, drag, similarity, open channels, pumps and networks.

  .. image:: _static/images/readme_fluid_engineering.png
     :alt: The Moody diagram, pump operating points, and open-channel specific energy curves
     :width: 100%

- :mod:`engineeringkit.heat_transfer` -- Conduction, convection, radiation, fins and heat exchangers.

  .. image:: _static/images/readme_heat_transfer.png
     :alt: Steady 2-D conduction in a plate, Planck blackbody spectra, and heat exchanger effectiveness-NTU curves
     :width: 100%

- :mod:`engineeringkit.materials` -- Crystals and diffusion, elastic constants, hardness, yield criteria, composites and laminates, fatigue and notches, fracture, creep and dislocation plasticity.

  .. image:: _static/images/readme_materials.png
     :alt: Ramberg-Osgood stress-strain curves, Tresca and von Mises yield surfaces, and Paris-law fatigue crack growth
     :width: 100%

- :mod:`engineeringkit.mechanical_engineering` -- Rigid-body motion, Lagrange's equations, machine elements, mechanisms, cams and multi-degree-of-freedom vibration.

  .. image:: _static/images/readme_mechanical_engineering.png
     :alt: Vibration transmissibility, normal modes of a five-mass chain, and a four-bar coupler curve
     :width: 100%

- :mod:`engineeringkit.reliability` -- Lifetimes and life data, accelerated testing, networks, structural reliability, FMEA, fault trees and maintenance.

  .. image:: _static/images/readme_reliability.png
     :alt: The bathtub hazard curve, a Weibull probability plot, and a Kaplan-Meier survival estimate
     :width: 100%

- :mod:`engineeringkit.signal_processing` -- Sampling and rate conversion, the FFT and DCT, spectra and AR models, IIR and FIR design, adaptive filters.

  .. image:: _static/images/readme_signal_processing.png
     :alt: The Gibbs phenomenon, Butterworth, Chebyshev and elliptic low-pass filters, and an STFT spectrogram of a chirp
     :width: 100%

- :mod:`engineeringkit.structural_engineering` -- Members, sections, beams on any supports, columns, stress and strain at a point, pressure vessels, energy methods, trusses, frames, cables and arches.

  .. image:: _static/images/readme_structural_engineering.png
     :alt: Pratt truss member forces, a two-span continuous beam's moment and deflection, and portal frame sway
     :width: 100%

- :mod:`engineeringkit.thermodynamics` -- Ideal and real gases, mixtures, the two laws, steady-flow devices, engine and refrigeration cycles, exergy, steam properties, psychrometrics and combustion.

  .. image:: _static/images/readme_thermodynamics.png
     :alt: An ideal Rankine cycle on the steam dome, a dual-cycle p-v diagram, and air-standard cycle efficiencies
     :width: 100%

Start here
----------

:doc:`quickstart` explains installation and basic use, and the
`quickstart notebook <https://colab.research.google.com/github/cpoli/engineeringkit/blob/main/notebooks/quickstart.ipynb>`__
runs in Colab or `in JupyterLite <lite/lab/index.html?path=quickstart.ipynb>`__.
Four :doc:`tutorials </tutorials/index>` follow one problem across several
domains: :doc:`/tutorials/beam_to_fatigue_life`,
:doc:`/tutorials/feedback_time_and_frequency`,
:doc:`/tutorials/power_plant_to_cooling_tower` and
:doc:`/tutorials/radio_link_to_bit_errors`. Each domain also has a set of
:doc:`exercises </exercises/index>` with worked solutions. The domain guides describe
units, sign conventions and model assumptions; :doc:`validation` records the
benchmark policy, and :doc:`limitations` identifies remaining scope limits.
Version 0.1.0 includes all twelve domains and enforces 100% production statement
coverage. That metric does not prove scientific correctness or model applicability.

.. code-block:: python

   from engineeringkit.structural_engineering import SimplySupportedBeam
   beam = SimplySupportedBeam(5, 200e9, 8e-6)
   result = beam.solve_point_load(10000)
   print(result.summary())
   fig, axes = result.plot()

.. toctree::
   :maxdepth: 2
   :caption: Subpackages
   :hidden:

   subpackages/index

.. toctree::
   :maxdepth: 2
   :caption: History
   :hidden:

   history/index

.. toctree::
   :maxdepth: 2
   :caption: Examples
   :hidden:

   examples/index

.. toctree::
   :maxdepth: 2
   :caption: API
   :hidden:

   api/index

.. rubric:: Related toolkits

Sibling Python toolkits for learning and teaching, built the same way:

- `blockchainkit <https://cpoli.github.io/blockchainkit/>`__ -- cryptography, consensus and blockchain systems
- `chemistrykit <https://cpoli.github.io/chemistrykit/>`__ -- computational chemistry
- `financialkit <https://cpoli.github.io/financialkit/>`__ -- quantitative finance
- `mathematicskit <https://cpoli.github.io/mathematicskit/>`__ -- computational mathematics
- `physicskit <https://cpoli.github.io/physicskit/>`__ -- computational physics
- `tbkit <https://cpoli.github.io/tbkit/>`__ -- tight-binding models, band topology and quantum transport
