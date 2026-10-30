Examples
========

.. include:: /_generated/nav/heat_transfer.rst

This gallery walks through every public feature of
``engineeringkit.heat_transfer``, one breakthrough per script: Newton's
cooling and the convection correlations, transient conduction from Biot to
Heisler, Fourier's law, fins and finite differences, blackbody radiation and
view factors and the radiation network, phase change, and heat exchangers.

:doc:`Units, assumptions and numerical validation </guides/heat_transfer>`
gives the physical validity limits. Each script is self-contained and can be
run directly with ``python examples/heat_transfer/<section>/<script>.py``.

Sections
--------

- **convection** -- Newton's law of cooling, the Dittus-Boelter and Churchill-Chu correlations, Graetz's tube flow, Pohlhausen's flat plate, the Chilton-Colburn analogy, Žukauskas's tube banks and Churchill and Bernstein's cylinder.
- **transient** -- Biot's number, Kelvin's age of the Earth and Heisler's charts.
- **conduction** -- Fourier's law, Harper and Brown's fins and the Crank-Nicolson scheme.
- **radiation** -- Kirchhoff's law, Stefan-Boltzmann, Wien's displacement, Nusselt's view factors and Oppenheim's network.
- **phase_change** -- Stefan's ice problem, Nusselt's film condensation, the boiling curve and Rohsenow's nucleate boiling.
- **exchangers** -- Kays and London's effectiveness-NTU method and Bowman, Mueller and Nagle's correction factor.
