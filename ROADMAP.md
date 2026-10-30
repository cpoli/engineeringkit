# Roadmap

What engineeringkit does not do yet, as a list of independent work items.
engineeringkit is a teaching toolkit, so the measure is not feature parity
with python-control, CoolProp or a finite element code. It is whether each
domain covers the canonical topics of an undergraduate course in that
field, tells their history, and checks that every example still shows the
right engineering.

## Where engineeringkit stands

Checked on 2026-10-20 against commit `854575c` plus the uncommitted working
tree.

- 0.1.0 is on PyPI (2026-10-18), with a GitHub Release, docs at
  https://cpoli.github.io/engineeringkit/ and the version DOI
  10.5281/zenodo.23452015 (concept record 23452014).
- 725 tests pass in about 5 s. Production **statement** coverage is 100%
  and enforced. **Branch** coverage is 99% (936 branches, 44 partial) and
  not measured in CI.
- `python -m mypy` is clean on the whole package (341 files), and the CI
  `typecheck` job runs exactly that.
- Twelve domains, all rebuilt to the `core/ systems/ utils/ visualizers/
  tests/` layout: 213 dated breakthroughs, one gallery example each.
- No domain has been checked against a course syllabus. The breakthrough
  lists were chosen for history, not for coverage, so canonical topics are
  missing in every domain (section D).

| Domain | Breakthroughs | Examples | Exercises |
|---|---:|---:|---:|
| `aerodynamics` | 27 | 29 | 3 |
| `communications` | 30 | 32 | 3 |
| `control_systems` | 26 | 26 | 3 |
| `electrical_engineering` | 24 | 24 | 3 |
| `fluid_engineering` | 25 | 25 | 3 |
| `heat_transfer` | 25 | 25 | 3 |
| `materials` | 26 | 26 | 3 |
| `mechanical_engineering` | 24 | 25 | 3 |
| `reliability` | 26 | 26 | 3 |
| `signal_processing` | 28 | 28 | 3 |
| `structural_engineering` | 29 | 29 | 3 |
| `thermodynamics` | 25 | 32 | 3 |

The history pages in `docs/source/history/` are the record of every entry.
The per-domain rebuild plans that used to be here are in git history
(`git show 854575c:ROADMAP.md`).

## How to work from this file

Each item below is meant to be done on its own, in one session, from a
prompt such as "Do ROADMAP item D7". An item names what to change and when
it is done. Where one item should come after another it says **After**.
Nothing else is ordered, but section A is small and should go first.

Every item is done when:

- `ruff check .`, `ruff format --check .`, `python -m mypy`,
  `MPLBACKEND=Agg pytest --cov=engineeringkit` and `make -C docs html`
  (which is `-W --keep-going`) are green;
- the change has a line under `## [Unreleased]` in CHANGELOG.md;
- the item is moved from "Open work" to "Done" here, with the date.

The project conventions apply to everything: SI inputs validated in
`_validation.py`, NumPy-style docstrings with units and the guide Notes
block, an analytical or benchmark test for every new model, 100%
production statement coverage, `filterwarnings = ["error"]`, the public
API kept backward compatible through each domain's `__init__.py`, and,
for each new breakthrough, a dated history entry with its own gallery
example, titled after it, that ends in a check of the result it
illustrates and is linked from no other entry. A history page must never
fall below 15 entries.

Sessions may run in parallel. CHANGELOG.md, this file, the table above,
`docs/source/conf.py` and `docs/source/validation.rst` are shared: re-read
them just before editing and change only your item's lines.

## Open work

### A. Repository hygiene

**A1. Commit or revert the pending working tree.** `docs.yml` and
`release.yml` now run `cd docs && make html` and publish `docs/build/html`,
but `docs/Makefile` and `docs/make.bat` are untracked. Pushing the
workflows without them breaks the docs deploy and the release job. Review
the diff (it also changes CLAUDE.md, AGENTS.md, README.md, VERIFICATION.md,
`docs/README.md`, `conf.py`, `quickstart.rst`, `pyproject.toml` and four
aerodynamics examples), commit what is intended, and drop the rest.
*Done when* `git status` is clean and the Docs workflow is green on `main`.

**A2. Stop generated plot images landing in the source tree.** Building the
docs writes `<page>-N.png` files from `.. plot::` blocks into
`docs/source/exercises/` and `docs/source/tutorials/`. Find out why
(`plot_html_show_source_link`, `plot_working_directory` or the build
target) and stop it; if that is not possible, add them to `.gitignore`.
*Done when* a clean docs build leaves `git status` clean.

**A4. Remove stale release and verification text.** `limitations.rst`
("Release status": nothing on PyPI, Colab and JupyterLite links do not
work), `VERIFICATION.md` (398 tests, 1878 statements, "JupyterLite launch
is intentionally not configured for this unpublished package") and any
README or guide text written before 0.1.0 describe a state that no longer
holds. Rewrite `VERIFICATION.md` as a short reproducible check list with
today's numbers, or fold it into `docs/source/validation.rst` and delete
it. *Done when* `grep -ri "unpublished\|not been released\|first release"`
finds nothing stale.

**A5. Add an `## [Unreleased]` section to CHANGELOG.md** for the two
commits after `v0.1.0` (dated example titles; Linux-only exercise check),
so the release workflow and the "Done" rule above have somewhere to write.

### B. Quality gates

**B1. Branch coverage.** Set `branch = true` under `[tool.coverage.run]`,
cover the 44 partial branches with real tests (no pragmas), and keep
`fail_under = 100`. A branch that cannot be reached is a sign to simplify
the code. *Done when* CI enforces 100% branch coverage.

**B2. Test the supported dependency range.** `pyproject.toml` promises
NumPy ≥ 1.24, SciPy ≥ 1.10 and Matplotlib ≥ 3.6, but CI only installs the
latest. Add a job on Python 3.10 with those minimum versions pinned, and
a `windows-latest` row to the test matrix (students use Windows). Fix
what fails or raise the floors. *Done when* both jobs are green.

**B4. Textbook reference problems, one domain per session.** 0.1.0 named
"reference problems" as the next step after statement coverage, and it
was never scheduled. For the domain named in the prompt, add
`engineeringkit/<domain>/tests/test_reference_problems.py` with at least
five worked examples from the textbooks listed for that domain in
section D, each citing author, edition and example number, and checking
the published answer to the precision the book prints. A disagreement is
either a bug to fix or a difference of method to document in the guide.
Twelve sub-items: B4.aerodynamics ... B4.thermodynamics.

**B5. Warn outside a model's range of validity.** `limitations.rst` says
assumptions "are stated, not enforced": Dittus-Boelter at Re = 500,
Colebrook in laminar flow or Prandtl-Glauert at M = 0.95 return a number
without a word. Add a `ValidityWarning(UserWarning)` to `_validation.py`
and a helper that emits it when an input leaves the range the docstring
states. Pilot it in `heat_transfer` and `fluid_engineering` (correlations
first), with `pytest.warns` tests, and document in CONTRIBUTING.md how to
silence it. New warnings are a behaviour change: record it in the
CHANGELOG. *Done when* the pilot domains warn and the other ten each have
a follow-up line (B5.<domain>) under this item.

### C. API consistency (before 1.0)

**C1. Deprecation policy and helper.** There is none, so no name can be
changed safely. Write the policy in CONTRIBUTING.md (deprecated names keep
working for at least one minor release and emit `DeprecationWarning`) and
add a small `_deprecated(old, new)` helper with a test. **Before** C2.

**C2. Naming audit.** List public names that break the pattern the rest of
the package follows, for example `rule_of_mixtures_E_longitudinal`
(capital), `db` (communications) next to `db20` (signal processing),
`natural_frequency` and `natural_frequency_hz`. Write the naming rule in
CONTRIBUTING.md, add the new names, and deprecate the old ones with C1's
helper. **After** C1.

**C3. One implementation per concept.** Several names exist in two
domains: `reynolds_number`, `dynamic_pressure` and `mach_number`
(aerodynamics, fluid_engineering), `hooke_stress` and `thermal_strain`
(materials, structural_engineering), `parallel_axis` (mechanical,
structural), and constants such as `G`/`STANDARD_GRAVITY`, `R_AIR`,
`BOLTZMANN` and `SPEED_OF_LIGHT` next to `engineeringkit.constants`.
Where two copies differ in signature or value, decide which is right;
keep one implementation and re-export it, and add a test asserting the
re-exports are the same object.

**C5. Result-object contract test.** CLAUDE.md says every result dataclass
is frozen and has `.summary()` (a dict of scalars) and `.plot()` returning
`(fig, ax)` or `(fig, axs)`, but no cross-domain test checks it. Add one
that discovers every `*Result` class through the domain `__all__` lists.

### D. Curriculum gaps, one domain per session

For the domain named in the prompt: compare its API and guide against the
table of contents of the two textbooks named, add the models listed below
(and any other canonical topic the comparison finds), each with an
analytical or benchmark test and a guide section. Where a model carries a
datable breakthrough, add a history entry and gallery example; check
dates and attributions against the primary source while writing the
entry, and say so in the entry when only a textbook could be found. Update
the table above. The lists were checked against the domain `__all__` on
2026-10-20; nothing in them exists yet.

All twelve items are done; see Done.

### E. History and citations

**E1. Citation audit.** The rebuild plans asked for every date and
attribution to be checked against the primary source, but nothing records
which ones were. Write a script (kept in `docs/` or the tests) that
extracts each entry's reference from the twelve history pages and flags
those with no DOI, stable URL or page range. Fix what can be fixed, and
list the rest in `limitations.rst` next to Jacobi and Thoma. Try once
more for those two: Jacobi's 1840s papers in *Annalen der Physik* or the
St Petersburg *Bulletin*, and Thoma's 1924 paper in the *Zeitschrift des
VDI* or the Munich hydraulics reports.

**E2. Test the exercise "Builds on" links.** Each link in
`docs/source/exercises/<domain>.rst` should point at a specific gallery
example, not a domain gallery index, and that example must exist. Add the
check to `engineeringkit/tests/test_history_conventions.py`.

### F. Course material

**F1. A tutorial for `aerodynamics`.** No tutorial uses it. Candidate:
*an aircraft from atmosphere to range* (`aerodynamics`, `thermodynamics`
for the Brayton engine, `structural_engineering` for the wing spar).

**F2. A tutorial for `electrical_engineering` and
`mechanical_engineering`.** Neither is used by a tutorial. Candidate:
*a motor drive from supply to shaft* (`electrical_engineering`,
`mechanical_engineering`, `control_systems`).

**F3. More exercises, one domain per session.** Three problems per domain
is a week of a course, not a term. Bring the named domain to eight,
ordered by difficulty, with at least one open-ended design problem, all
checked by `test_course_material.py` as today.

**F4. Syllabus maps.** Instructors need to see what the toolkit covers in
the course they teach. Add `docs/source/syllabus/<domain>.rst` mapping
each chapter of the two textbooks in section D to the functions, gallery
examples and exercises that cover it, with the uncovered chapters left
visible. It can be written before section D and updated by it.

### G. Releases

**G1. Release 0.2.0.** **After** section A. Bump the version in
`engineeringkit/__init__.py`, `pyproject.toml` and `CITATION.cff`, rename
`## [Unreleased]` in CHANGELOG.md, tag `v0.2.0`, then add the release's
own Zenodo version DOI (not the concept DOI) to `CITATION.cff` and the
README badge.

**G2. conda-forge.** Submit a `recipe.yaml` to `conda-forge/staged-recipes`
built from the PyPI sdist, with only `version` in `context`, and
`name: engineeringkit` and the source URL written out. Add the conda
badge and install line to the README when it is merged.

### H. Scope decisions

- **Out of scope**, since the dependencies stay NumPy, SciPy and
  Matplotlib: CAD and meshing, general finite element or CFD solvers,
  wrappers around CoolProp, python-control or SPICE, and unit-carrying
  types (inputs stay raw SI numbers).
- **Sibling kits.** physicskit has `fluids` and `classical`, chemistrykit
  has `thermo` and `kinetics`. Nothing records where engineeringkit stops
  and they start. **H1.** Compare the overlapping APIs, write the boundary
  in `docs/source/limitations.rst` (engineering correlations and design
  here; first-principles physics and chemistry there), and cross-link the
  matching gallery examples.
- **Candidate domains**, not scheduled, each needing 15+ breakthroughs and
  several sessions: `geotechnical` (Terzaghi's effective stress and
  consolidation, Mohr-Coulomb, Rankine earth pressure, Casagrande),
  `power_systems` (see D4), `manufacturing` (Taylor tool life, Merchant's
  cutting force, statistical process control). Start one only after
  section D, since every existing domain still has gaps.

## Done

Finished items move here with the release they shipped in, as in
CHANGELOG.md.

### Unreleased

- **B3. External links checked weekly** (2026-10-25). The
  `Linkcheck` workflow runs `make -C docs linkcheck` every Monday and on
  demand, without executing the gallery. The 14 history-page DOIs are
  now `:doi:` links (`sphinx.ext.extlinks`), so they are checked too.
  `linkcheck_ignore` skips the JupyterLite site and, by DOI prefix, the
  publishers that answer bots with 403: Wiley, ACS, ASCE, the Royal
  Society and ACM. Local run: 16 links working, 0 broken; the first
  scheduled run is still to come.
- **A3. CI type check covers the whole package** (2026-10-25). The
  `typecheck` job runs `python -m mypy` with the `pyproject.toml`
  settings instead of `mypy engineeringkit/integrators`, and checks all
  341 files (282 when the item was written).
- **D1. `aerodynamics` curriculum gaps** (2026-10-23). Against
  Anderson's *Fundamentals of Aerodynamics* and *Introduction to Flight*:
  oblique shocks with Thompson's closed-form wave angle, weak and strong
  branches and maximum deflection (`ObliqueShockResult`); the
  Prandtl-Meyer function, its inverse and expansion fans; Fanno and
  Rayleigh duct flow with the choking length and heat (`DuctFlowState`,
  Darcy friction factor); a general lifting-line solver by Glauert's
  Fourier series for any chord, twist and zero-lift distribution
  (`LiftingLineResult`); the Hess-Smith source-vortex panel method
  (`PanelResult`); the thin-airfoil quarter-chord moment and the
  aerodynamic centre; tail volume, downwash gradient, neutral point,
  moment slope and static margin; rate and angle of climb, load factor,
  turn radius and rate, take-off ground roll, and propeller Breguet
  range. Six new breakthroughs (22 to 27) and second examples on the
  thin-airfoil and Breguet entries. The Fanno thesis and Meyer's and
  Prandtl's papers could not be checked; those entries cite Anderson's
  *Modern Compressible Flow* and say so. Still missing from those
  textbooks and left for later: shock reflection and interaction,
  conical shocks (Taylor-Maccoll), method of characteristics and nozzle
  design, unsteady waves and shock tubes; vortex and linear-strength
  panel methods, viscous-inviscid coupling and airfoil drag, swept and
  delta wings, lifting-surface and vortex-lattice methods; transonic
  critical Mach number and drag divergence; dynamic stability (phugoid,
  short period), trim and control, and elevator and stick-free stability;
  endurance, landing distance, energy-height methods and V-n diagrams;
  and the hypersonic Newtonian theory.
- **D6. `heat_transfer` curriculum gaps** (2026-10-23). Against
  Incropera et al. and Çengel and Ghajar: the flat plate in laminar,
  turbulent and mixed flow (Pohlhausen); Churchill-Bernstein for the
  cylinder and Žukauskas for tube banks with the row correction; laminar
  tube Nusselt numbers, entry lengths and Hausen's entry-region average;
  the Chilton-Colburn heat and mass transfer analogy; fins with a
  convecting tip; steady 2-D conduction by finite differences
  (`Conduction2DResult`) and four conduction shape factors; Heisler's
  cylinder and sphere series; Oppenheim's gray-surface network
  (`EnclosureResult`) with reradiating surfaces, the two-surface form,
  radiation shields and view-factor reciprocity and summation; crossflow
  and shell-and-tube effectiveness and Bowman, Mueller and Nagle's LMTD
  correction factor; Rohsenow's nucleate boiling. Eight new breakthroughs
  (18 to 25): Graetz 1883, Pohlhausen 1921, Chilton-Colburn 1934,
  Bowman-Mueller-Nagle 1940, Rohsenow 1952, Oppenheim 1956, Žukauskas
  1972, Churchill-Bernstein 1977. Hausen (1943) and the tube-bank
  constants could be checked only through Incropera et al., and the
  entries say so. Still missing from those textbooks and left for later:
  conduction with heat generation and contact resistance, multidimensional
  transients by product solutions, annular fins and fin arrays, turbulent
  internal flow beyond Dittus-Boelter (Gnielinski appears only in an
  example), non-circular ducts and annuli, free convection on horizontal
  plates, cylinders, spheres and in enclosures, jets and packed beds, film
  boiling, flow boiling and condensation on tubes, spectral and solar
  radiation (band fractions, emissivity data) and gas radiation, compact
  exchangers and fouling factors, and diffusion mass transfer beyond the
  analogy.

- **D8. `mechanical_engineering` curriculum gaps** (2026-10-23). Against
  Shigley and Rao: normal modes from the M-K eigenproblem (`ModalResult`),
  closed-form modal superposition with modal damping (`MDOFResult`) and
  the direct harmonic response; the forced SDOF magnification factor and
  phase; Rayleigh's quotient and Dunkerley's sum; Lewis gear-tooth bending
  with the tabulated form factor and Barth's velocity factor; solid shafts
  in combined bending and torsion (von Mises, Tresca, DE-Goodman); helical
  springs (rate, index, Wahl and Bergsträsser factors); bolted joints
  (bolt and frustum member stiffness, `BoltedJoint` with load factor and
  joint diagram, preload torque); one- and two-plane balancing
  (`BalanceResult`); four-bar and slider-crank velocity and acceleration
  (`FourBarMotion`, `SliderCrankMotion`); six cam displacement laws
  (`CamMotion`); Lagrange's equations by numerical differentiation of a
  user Lagrangian. Eight new breakthroughs (16 to 24) and a second
  Rayleigh example on forced vibration. Still missing from those
  textbooks and left for later: base excitation, Duhamel's integral and
  shock spectra, general (non-proportional) damping, continuous systems
  (strings, bars, Euler-Bernoulli beam modes), Newmark integration,
  vibration measurement and nonlinear and random vibration (Rao); AGMA
  gear stresses, helical, bevel and worm gears, planetary trains, welded
  joints, clutches and brakes, chains, spring buckling and surge, and
  cam profile synthesis with pressure angle (Shigley and kinematics
  texts).

- **D4. `electrical_engineering` curriculum gaps** (2026-10-23). Scope
  decision first: the load flow, swing equation and fault models stay in
  this domain, since section H defers new domains until section D is done
  and a `power_systems` domain would need 15 breakthroughs of its own.
  Against Hayt, Kemmerly and Durbin and Glover, Sarma and Overbye: mesh
  analysis (`MeshResult`); the RL current and the over-, critically and
  underdamped series RLC step; ideal inverting, non-inverting and
  integrating op amps; the Shockley diode and its load line; DC machine
  and induction motor torque-speed curves with breakdown torque; the
  per-unit system; fault currents for all four fault types from the
  sequence impedances; Newton-Raphson load flow with PV buses; the swing
  equation (`SwingResult`), equal-area critical clearing angle and
  critical clearing time. Eight new breakthroughs (17 to 24); the 1926
  Evans and Wagner stability paper could not be checked and the entry
  says so. Still missing from those textbooks and left for later:
  supermeshes and dependent sources, AC (phasor) nodal and mesh analysis,
  coupled inductors and the non-ideal transformer equivalent circuit,
  two-port parameters, Laplace-domain circuit analysis and circuit
  frequency response, Y-delta conversion, op-amp non-idealities,
  transistor circuits, ABCD long-line models, the synchronous machine,
  the bus impedance matrix and fault studies built on it, economic
  dispatch, PV-to-PQ switching at reactive limits and multimachine
  stability.
- **D10. `signal_processing` curriculum gaps** (2026-10-23). Against
  Oppenheim and Schafer and Proakis and Manolakis: Chebyshev type I and II
  and elliptic filters with minimum-order formulas; the bilinear transform
  written out, with prewarping; Parks-McClellan equiripple FIR design with
  Kaiser's length estimate; circular convolution, zero padding and
  overlap-add; downsampling, upsampling, decimation and interpolation; the
  Hilbert transform, analytic signal, envelope, instantaneous frequency and
  an FIR Hilbert transformer; the Levinson-Durbin recursion, Yule-Walker
  fitting and the AR spectrum (`ARModel`); the orthonormal DCT; the real
  cepstrum. Eleven new breakthroughs (17 to 28). Still missing from those
  textbooks and left for later: impulse invariance, high-pass, band-pass
  and band-stop transformations of the new designs (only low-pass), Bessel
  filters, lattice and direct-form II realizations and coefficient
  quantization, polyphase filters and rational resampling as one
  operation, overlap-save, the chirp z-transform, Burg's method and
  model-order criteria (AIC), the minimum-phase/all-pass decomposition and
  the complex cepstrum, and multirate filter banks.

- **D3. `control_systems` curriculum gaps** (2026-10-23). Against Ogata
  and Franklin, Powell and Emami-Naeini: lead, lag and lead-lag
  compensators with Ogata's frequency-response design procedures
  (`CompensatorResult`); system type and the static error constants;
  sensitivity, complementary sensitivity and the margins guaranteed by
  `Ms`; the Nichols chart with M-contours; Tustin's bilinear transform
  with prewarping; PID with a saturating actuator and clamping or
  back-calculation anti-windup; numerical and closed-form describing
  functions (relay with hysteresis, saturation, dead zone), limit-cycle
  prediction with Loeb's stability test, and the phase plane
  (`PhasePortraitResult`). Five new breakthroughs (21 to 26): Andronov,
  Vitt and Khaikin 1937, Bode's sensitivity integral 1945, the Nichols
  chart 1947, Tustin 1947, Kochenburger 1950. Still missing from those
  textbooks and left for later: root-locus design of compensators (angle
  and magnitude conditions), the Bode gain-phase relation as a function,
  PID derivative filtering and two-degree-of-freedom set-point weighting,
  discrete-time design (direct z-plane, deadbeat), MIMO frequency
  response and singular values, and Popov or circle criteria.
- **D7. `materials` curriculum gaps** (2026-10-23). Against Callister
  and Dowling: cubic lattice parameters, packing factors, theoretical
  density, interplanar spacing, Bragg angles and powder-pattern selection
  rules; Fick's first law, the Arrhenius diffusivity and the
  error-function carburizing profile and time; the lever rule and Avrami
  kinetics with a fit; Brinell and Vickers hardness and the
  hardness-strength conversion; Gerber and Soderberg safety factors next
  to Goodman, the endurance-limit estimate and Marin factors; Inglis's
  elliptical hole, Peterson and Neuber notch sensitivity, the fatigue
  notch factor and Neuber's rule on a Ramberg-Osgood curve;
  three-dimensional principal and von Mises stresses; Halpin-Tsai and
  classical laminate theory (`Laminate`); Weibull brittle strength built
  on `reliability.weibull_cdf`. Eight new breakthroughs (18 to 26). Still
  missing from those textbooks and left for later: phase diagrams as data
  (only the lever rule on a given tie line), TTT/CCT diagrams,
  strengthening by solid solution and precipitation, polymers, ceramics
  and corrosion, electrical, thermal and magnetic properties (Callister
  part 2), cyclic stress-strain curves, mean-stress strain-life
  (Morrow, Smith-Watson-Topper), the J-integral and CTOD, and laminate
  failure criteria (Tsai-Wu, Hashin).
- **D5. `fluid_engineering` curriculum gaps** (2026-10-23). Against
  White and Munson et al.: the momentum equation for jets, vanes, bends
  and Pelton's wheel; pump and system curves, the operating point in
  series and parallel, the affinity laws and Euler's turbomachine head;
  specific energy, critical and alternate depths, Bélanger's hydraulic
  jump; orifices, tank draining, rectangular and V-notch weirs with
  Rehbock's coefficient; Stokes drag, the sphere drag curve and terminal
  velocity; Couette flow and the Couette viscometer. Eight new
  breakthroughs (17 to 25). Still missing from those textbooks and left
  for later: boundary layers inside this domain (only `aerodynamics` has
  the laminar one), potential flow, compressible duct flow (in
  `aerodynamics`), turbine selection charts and gradually varied flow
  profiles.
- **D2. `communications` curriculum gaps** (2026-10-24). Against Haykin,
  *Communication Systems*, and Proakis and Salehi: North's matched filter
  and its 2E/N0 peak SNR; Gray-labelled M-PSK and square M-QAM with
  minimum-distance detection and exact symbol-error rates (Craig's
  integral for PSK); raised-cosine and root-raised-cosine pulses, their
  spectrum and bandwidth, and eye diagrams; PCM with its SQNR, mu-law
  companding and linear delta modulation with the slope-overload limit;
  general entropy and Huffman coding (`HuffmanCode`); CRCs by polynomial
  division with CRC-16-CCITT; convolutional encoding, free distance and
  hard- and soft-decision Viterbi decoding (`ViterbiDecodeResult`);
  Rayleigh fading, MRC and selection diversity; effective aperture,
  aperture gain, directivity from a pattern and the half-wave dipole; and
  the output SNR of AM and FM with the FM threshold. Thirteen new
  breakthroughs (18 to 30) and second examples on the Nyquist and Friis
  transmission entries. Reeves's exact 1938 filing date, the co-inventors
  of the 1946 delta-modulation patent, the page range of Carter's 1932
  paper and the first published 73-ohm calculation could be found only in
  secondary sources; the entries say so. Still missing from those
  textbooks and left for later: SSB and VSB generation and the
  superheterodyne receiver; FM pre-emphasis; the phase-locked loop and
  carrier and symbol-timing recovery; FSK, MSK, DPSK and noncoherent
  detection; orthogonal signalling and union bounds; equalization
  (zero-forcing, MMSE, decision feedback) and partial-response
  signalling; adaptive delta modulation, DPCM and sigma-delta; Lempel-Ziv
  coding and rate-distortion; BCH, Reed-Solomon, turbo and LDPC codes;
  Rician fading, delay spread, coherence bandwidth and Doppler (Clarke and
  Jakes); spread spectrum and CDMA; and MIMO.
- **D12. `thermodynamics` curriculum gaps** (2026-10-24). Against Çengel
  and Boles, *Thermodynamics: An Engineering Approach*, and Moran and
  Shapiro: steady-flow energy balances for nozzles, turbines,
  compressors, throttles and mixing chambers with isentropic
  efficiencies; `rankine_cycle` with superheat, reheat, an open
  feedwater heater and real turbines and pumps (`RankineResult`), which
  reproduces Çengel and Boles's three 15 MPa plants; `brayton_cycle` with
  real components, regeneration, intercooling and reheat
  (`BraytonResult`); the Stirling and dual cycles; the Peng-Robinson
  equation of state, its saturation pressure and R-134a properties
  (`CubicFluid`, `R134A`, `refrigerant_properties`), which close the
  refrigerant gap in `limitations.rst`, and `vapor_compression_cycle`
  (`VaporCompressionResult`); Dalton's ideal-gas mixtures and the entropy
  of mixing; the compressibility factor, Pitzer's acentric factor and the
  generalized virial correlation; Maxwell's relations as `cp_minus_cv`
  and the general Joule-Thomson coefficient; moist-air enthalpy,
  adiabatic saturation, mixing, cooling coils and cooling towers;
  combustion stoichiometry, Hess's-law heating values and the adiabatic
  flame temperature on NIST Shomate enthalpies; IAPWS-IF97 region 5
  (2007 revision, to 50 MPa) and an array-accepting `steam_properties`.
  Six new breakthroughs (20 to 25): Dalton 1801, Stirling 1816, Hess
  1840, Maxwell 1871, Pitzer 1955, Peng and Robinson 1976, and second
  examples on the Joule, Perkins, Rankine, Brayton, Diesel, Carrier and
  IF97 entries. Dalton's 1802 volume, Stirling's patent and Maxwell's
  1871 first edition could not be checked online; those entries say so,
  and Hess's citation is corrected to *Ann. Phys.* 126 with its DOI. Still
  missing from those textbooks and left for later: transient
  (unsteady-flow) tank filling and discharge; closed feedwater heaters,
  cogeneration and combined gas-steam cycles; jet-propulsion and Ericsson
  cycles as functions; cascade, multistage and gas refrigeration and
  absorption chillers; real-gas departure charts (Lee-Kesler) and
  reference-quality refrigerant equations; real-gas and liquid mixtures;
  evaporative coolers and heating with humidification as named
  processes; rich combustion, dissociation and chemical equilibrium
  (Gibbs minimization, equilibrium constants); fuel cells; and exergy
  balances of the new cycles.
- Gallery example titles carry their breakthrough's date.
- An exercise check that only passed on Linux is fixed.
- **C4. Array-input policy** (2026-10-23). Policy in CONTRIBUTING.md.
  `engineeringkit/tests/test_array_inputs.py` finds every scalar function
  from the parameter types in its docstring, and calls it with 1-D arrays
  built from scalar calls recorded from its own domain tests. Of the
  scalar functions that existed at `854575c`, 21 were vectorized, 7 have
  documented scalar-only parameters and the rest already broadcast (`naca4_zero_lift_angle`,
  `plane_wall_transient`, `neumann_stefan_front`, `crack_growth_life`,
  `jeffcott_response`, `tuned_mass_damper`, `kaiser_order`).
  `steam_properties` broadcasts since the D12 work. The 39 functions that
  in-progress D items had added but not yet vectorized are in the test's
  `PENDING` list (strict xfail), which each item clears.
- **D11. `structural_engineering` curriculum gaps** (2026-10-24). Against
  Hibbeler's *Mechanics of Materials* and *Structural Analysis* and Gere
  and Goodno: `Beam`, an Euler-Bernoulli beam on any pin, roller and fixed
  supports solved by Macaulay's singularity functions with the reactions
  as unknowns (point loads, couples, uniform and linearly varying loads;
  determinate and indeterminate alike), and `superpose_beam_results`;
  `composite_section` with I, T, channel, angle and box sections and
  transformed sections by modular ratio (`SectionProperties`, with
  principal axes, Q and cut width); the generalized flexure formula and
  neutral-axis angle; Jourawski's VQ/It, shear flow and the channel shear
  centre; Bredt's thin-walled torsion; thin-walled vessels and Lamé's
  thick cylinder; three-dimensional principal stresses; strain rosettes
  and plane-stress Hooke's law; effective-length factors, the secant
  formula and its inverse, and Johnson's parabola; the unit-load method
  for trusses and beams; Maney's slope-deflection with chord rotations
  (`SlopeDeflectionResult`); parabolic and catenary cables and the
  parabolic three-hinged arch (`ArchResult`). Twelve new breakthroughs
  (17 to 29): Hooke's inverted chain 1675, the catenary 1691, Young 1807,
  Lamé 1833, Jourawski 1856, Mohr's unit load 1874, Johnson 1893, Bredt
  1896, Maney 1915, Macaulay 1919, Maillart 1921, Simmons and Ruge 1938.
  Mohr's 1882 entry and example gain the three circles. Only secondary
  accounts could be found for Hooke, the catenary papers, Young, Johnson,
  Maney and the strain-gauge dates, and Bredt's and Lamé's page numbers
  come from bibliographies; the entries say so. Eight of the new
  functions (the five shape builders, `catenary_cable`,
  `secant_formula_load`, `thick_wall_cylinder_stresses`) are scalar-only
  and xfailed as pending vectorization in `test_array_inputs.py`. Still
  missing from those textbooks and left for later: variable EI and
  internal hinges in `Beam`, the moment-area and conjugate-beam methods,
  shear deformation (Timoshenko beams) and beams on elastic foundations;
  open-section (Saint-Venant) torsion constants and warping; rotating
  discs and shrink fits; sway frames in slope-deflection, approximate
  frame analysis (portal and cantilever methods) and Müller-Breslau
  influence lines for indeterminate structures; cables under point
  loads, two-hinged and fixed arches; plates and shells; and design-code
  checks (load combinations, LRFD resistance factors).
- **D9. `reliability` curriculum gaps** (2026-10-24). Against Ebeling,
  *An Introduction to Reliability and Maintainability Engineering*, and
  O'Connor and Kleyner, *Practical Reliability Engineering*: lognormal and
  gamma lifetimes; Johnson's exact median ranks, Benard's approximation and
  adjusted ranks for censored data; Weibull rank regression (on y or x) and
  Cohen's censored maximum likelihood (`WeibullFitResult` with the B-life
  and a probability plot); Arrhenius, Eyring and inverse-power-law
  acceleration factors and an Arrhenius fit; minimal path and cut sets of
  general block diagrams with exact reliability and Esary-Proschan bounds;
  cold standby with an imperfect switch; FMEA risk priority numbers,
  ranking (`FMEAResult`) and MIL-STD-1629A criticality; Monte Carlo
  failure probability with importance sampling at the design point
  (`MonteCarloResult`); gamma-Poisson updating of a failure rate
  (`FailureRatePosterior`) and a prior from a median and upper bound. Nine
  new breakthroughs (18 to 26). Johnson's 1951 paper and the text of
  MIL-P-1629 could not be found; those entries cite later sources and say
  so. Still missing from those textbooks and left for later: confidence
  bounds on Weibull parameters (Fisher matrix, likelihood ratio), the
  three-parameter Weibull, and life-stress models fitted by maximum
  likelihood; warm standby, load sharing and common-cause (beta-factor)
  models; reliability allocation and test planning (sequential and
  Bayesian demonstration tests); spares provisioning and maintainability
  demonstration; importance measures for general structures; and
  second-order (SORM) and system-level structural reliability.

### 0.1.0 (2026-10-18)

- **Release infrastructure.** Public GitHub repository, `CITATION.cff`,
  tag-driven `release.yml` (version check against `__version__`,
  `pyproject.toml` and `CITATION.cff`, PyPI trusted publishing, GitHub
  Release from the CHANGELOG, `gh-pages` rebuild), `docs.yml` deploying
  `gh-pages` from `main`, CI over Python 3.10-3.14 on Linux and macOS,
  the Zenodo version DOI.
- **Blocking type check** (2026-10-13). `python -m mypy` checks the whole
  package with `check_untyped_defs`, `strict_equality` and
  `warn_unreachable`. CI ran it on `integrators` only until A3.
- **All twelve domains rebuilt** to the `core/ systems/ utils/
  visualizers/ tests/` layout, each with a dated history page of 16 to 21
  breakthroughs and one gallery example per breakthrough; the per-domain
  `plot_<domain>_models.py` and benchmark examples, `test_phase1.py` and
  the single-module layout descriptions were removed afterwards.
- **Numerical gaps** found in 0.1.0's `limitations.rst`: high-order IIR
  filters as second-order sections, stability margins over every
  crossover, truss and frame conditioning diagnostics, real fluids (van
  der Waals, Clausius-Clapeyron, Antoine, IAPWS-IF97 regions 1 to 4),
  complete NumPy-style docstrings, array `gamma` in the thermodynamics
  relations, and Castigliano integration along a beam.
- **Course material** (2026-10-13). Four tutorials across two or three
  domains each, an exercise set of three problems per domain, and
  `notebooks/quickstart.ipynb` for Colab and JupyterLite, all re-run and
  number-checked by `test_course_material.py`.
- Physical input contracts, independent numerical benchmarks in every
  domain, twelve assumptions/validity guides, and 100% production
  statement coverage enforced in CI.
