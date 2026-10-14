Tutorials
=========

The galleries show one idea at a time. Each tutorial follows one
engineering problem through several domains, the way a design or a
course project does. Every code block runs when the documentation is
built, and the test suite checks that the numbers printed under each one
are still what the code produces.

.. grid:: 1 2 2 2
   :gutter: 2

   .. grid-item-card:: A beam from load to fatigue life
      :link: beam_to_fatigue_life
      :link-type: doc

      Structural engineering, materials and reliability: moments and
      deflections, Goodman, Basquin and Miner, and a Weibull life for a
      four-beam runway.

   .. grid-item-card:: A feedback loop in the time and frequency domains
      :link: feedback_time_and_frequency
      :link-type: doc

      Control systems and signal processing: Ziegler-Nichols tuning, step
      response and phase margin, a frequency response measured from
      samples, and the phase cost of a sensor filter.

   .. grid-item-card:: A power plant from cycle to cooling tower
      :link: power_plant_to_cooling_tower
      :link-type: doc

      Thermodynamics, heat transfer and fluid engineering: a Rankine cycle,
      exergy, a condenser by effectiveness-NTU, and the cooling-water pump.

   .. grid-item-card:: A radio link from antenna to bit errors
      :link: radio_link_to_bit_errors
      :link-type: doc

      Communications and signal processing: a link budget, cascaded noise,
      a simulated matched filter, and a Hamming code.

To try the package without installing anything, open the
`quickstart notebook <https://colab.research.google.com/github/cpoli/engineeringkit/blob/main/notebooks/quickstart.ipynb>`__
in Colab or `in JupyterLite <../lite/lab/index.html?path=quickstart.ipynb>`__.
Both install engineeringkit from PyPI.

.. toctree::
   :hidden:

   beam_to_fatigue_life
   feedback_time_and_frequency
   power_plant_to_cooling_tower
   radio_link_to_bit_errors
