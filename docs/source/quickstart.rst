Installation and quickstart
===========================

Install from PyPI::

    python -m pip install engineeringkit

For development::

    python -m pip install -e '.[dev]'
    MPLBACKEND=Agg python -m pytest
    python -m ruff check .
    make -C docs html
    python -m build

Beam response
-------------

.. code-block:: python

    from engineeringkit.structural_engineering import SimplySupportedBeam
    result = SimplySupportedBeam(5, 200e9, 8e-6).solve_point_load(10000)
    print(result.summary())
    fig, axes = result.plot()

Damped vibration
----------------

.. code-block:: python

    import numpy as np
    from engineeringkit.mechanical_engineering import simulate_sdof
    result = simulate_sdof(1, .2, 4, x0=.1, t=np.linspace(0, 10, 501))
    print(result.summary())
    fig, ax = result.plot()
    fig, animation = result.animate()
    animation.save("oscillator.gif", writer="pillow", fps=30)

Examples
--------

Run from the repository after installation::

    MPLBACKEND=Agg python examples/numerical_methods.py
