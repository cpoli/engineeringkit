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

- 0.2.0 is on PyPI (2026-10-30), with a GitHub Release, docs at
  https://cpoli.github.io/engineeringkit/ and the version DOI
  10.5281/zenodo.23795071 (concept record 23452014). 0.1.0
  (2026-10-18) is 10.5281/zenodo.23452015.
- 1730 tests pass in about 10 s. Production statement and branch coverage
  (1856 branches) are 100% and enforced in CI.
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
production statement and branch coverage, `filterwarnings = ["error"]`,
the public API kept backward compatible through each domain's `__init__.py`, and,
for each new breakthrough, a dated history entry with its own gallery
example, titled after it, that ends in a check of the result it
illustrates and is linked from no other entry. A history page must never
fall below 15 entries.

Sessions may run in parallel. CHANGELOG.md, this file, the table above,
`docs/source/conf.py` and `docs/source/validation.rst` are shared: re-read
them just before editing and change only your item's lines.

## Open work

### A. Repository hygiene

All items are done; see Done.

### B. Quality gates

**B4. Textbook reference problems, one domain per session.** 0.1.0 named
"reference problems" as the next step after statement coverage, and it
was never scheduled. For the domain named in the prompt, add
`engineeringkit/<domain>/tests/test_reference_problems.py` with at least
five worked examples from the textbooks listed for that domain in
section D, each citing author, edition and example number, and checking
the published answer to the precision the book prints. A disagreement is
either a bug to fix or a difference of method to document in the guide.
Twelve sub-items: B4.aerodynamics ... B4.thermodynamics, all done (see
Done).

**B5. Warn outside a model's range of validity, one domain per session.**
The pilot is done (see Done): `_validation.validity()` issues an
`engineeringkit.ValidityWarning` when an input leaves the range a
docstring states, as CONTRIBUTING.md ("Range of validity") describes. For
the domain named in the prompt, find the models whose docstrings or guide
state a range that is not checked, wire `validity()` into those whose
inputs are enough to test it, add `pytest.warns` tests on both sides of
each range, fix the tests and gallery examples that now warn, and record
the new warnings in the CHANGELOG. A docstring search on 2026-11-01 found
these candidates:

- **B5.aerodynamics**: `prandtl_glauert_cp` above M = 0.7 (the docstring
  says unreliable there), `skin_friction_laminar` above Re = 5e5, and the
  lifting-line wing's aspect-ratio and sweep limits.
- **B5.communications**: the envelope-detector noise model at low
  carrier-to-noise ratio and FM output SNR below `fm_threshold_snr`.
- **B5.control_systems**: no stated range found by the search; audit the
  guide's assumptions (linearisation, Padé delay orders).
- **B5.electrical_engineering**: no stated range found by the search;
  audit short- and medium-line models against line length.
- **B5.materials**: `brinell_hardness` outside 0.24D < d < 0.6D, and
  linear-elastic fracture where the plastic zone is not small.
- **B5.mechanical_engineering**: first-cycle yielding in the shaft fatigue
  models (`shafts.py`).
- **B5.reliability**: no stated range found by the search; audit the
  acceleration models (Arrhenius, Peck) against their fitted ranges.
- **B5.signal_processing**: no stated range found by the search; audit the
  window and filter-design approximations.
- **B5.structural_engineering**: Euler buckling above the proportional
  limit (`stability.py`).
- **B5.thermodynamics**: no stated range found by the search; audit the
  ideal-gas models at high reduced pressure and the psychrometric
  temperature range.

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

- **B4.aerodynamics. Textbook reference problems** (2026-11-02).
  `aerodynamics/tests/test_reference_problems.py` checks ten worked
  examples from Anderson: *Fundamentals of Aerodynamics*, 3rd ed.,
  Examples 4.4 (thin-airfoil lift), 4.8 (turbulent skin friction), 5.2
  (induced drag at a new aspect ratio), 8.5 (stagnation and sonic
  state), 9.2 (oblique shock) and 9.8 (Prandtl-Meyer expansion), and
  *Modern Compressible Flow*, 3rd ed., Examples 3.4 (normal shock), 3.10
  (thermal choking), 3.11 and 3.13 (Fanno pipe flow and choking length).
  Example 3.13 found a bug: `fanno_flow`'s friction parameter and
  `fanno_choking_length` were four times too large, since the Darcy form
  applied the Fanning factor of four twice; the ODE benchmark and the
  gallery check had been written to agree with it, and are fixed too.
  Examples 9.2 and 9.8 differ only because Anderson rounds to the nearest
  table entry (M_n1 = 1.6, M2 = 2.0); the guide says so. The printed
  answers were taken from the FOSSEE Python Textbook Companion
  transcriptions of both books, not from the books themselves.

- **B4.reliability. Reference problems** (2026-11-02).
  `reliability/tests/test_reference_problems.py` checks eight worked
  examples of the NIST/SEMATECH e-Handbook of Statistical Methods, ch. 8,
  against the printed answers, since no legitimate copy of Ebeling or
  O'Connor and Kleyner was at hand: 8.1.5.1 (Arrhenius factors 133 and
  17,597), 8.4.5.1 (800 h, 2 failures: MTBF bounds 127 and 2251 h; 267 h
  with none), 8.4.1.3 (censored Weibull MLE, eta = 606.5280 h, beta =
  1.7208, AIC 154.27), 8.2.2.1 (probability plot, slope 1.46), 8.2.5
  (gamma prior a = 2.863, b = 1522.46 h), 8.4.6 (posterior MTBF bounds
  495, 600, 901 and 1897 h), 8.2.1.5 (Kaplan-Meier products) and 8.4.3
  (lognormal F = 0.014). All agree; nothing needed fixing. The guide
  records where the handbook's rounding moves a printed number: 17,597
  needs 1/k = 11605 K/eV and 273.16 K, the plot's 657 h comes from a
  rounded slope (the fit gives 667 h), and its MLE scale and prior rate
  stop short of the optima. Textbook examples from Ebeling and O'Connor
  and Kleyner can be added later.

- **B4.communications. Textbook reference problems** (2026-11-02).
  `communications/tests/test_reference_problems.py` checks six worked
  examples of Haykin, *Digital Communications* (Wiley, 1988), against
  the printed answers: Example 2.2 (H = 1.5 bits, 3 bits for the
  second-order extension), 2.3 (Huffman code, L = 2.2, H = 2.12,
  variance 0.16), 5.3 (8-bit PCM of 4 kHz, B = 32 kHz, 40.8 dB), 8.3
  (cyclic (7,4) code, g = 1 + X + X^3, message 1001 to 0111001), 8.5
  (syndromes 000 and 110) and 8.7 (convolutional encoder, 10011 to
  11 10 11 11 01 01 11). Section D names Haykin's *Communication
  Systems* and Proakis and Salehi, but neither text could be consulted,
  so the examples come from *Digital Communications*, whose example
  numbers, pages and answers were checked against the FOSSEE Python
  Textbook Companion (https://tbc-python.fossee.in/book-details/977/).
  All agree, so nothing was fixed. The guide records the one difference
  of method: Haykin's 6n - 7.2 dB loads the quantizer to four times the
  rms signal, 9.03 dB below `pcm_sqnr_db`'s full-scale sine, and rounds
  to 6 dB per bit (40.8 against 40.89 dB). Examples from the two
  section D books are left for whoever has copies.

- **B4.fluid_engineering. Textbook reference problems** (2026-11-02).
  `fluid_engineering/tests/test_reference_problems.py` checks five worked
  examples against the printed answers. White, *Fluid Mechanics*, 5th
  ed.: Example 6.7 (oil in a sloping cast-iron pipe, V = 6.37 m/s,
  Re = 127,000, h_f = 117 m), 6.9 (V = 4.84 m/s, Q = 0.342 m^3/s from an
  8 m head loss), 6.10 (the same worked back to d = 0.300 m) and 10.8
  (hydraulic jump, Fr1 = 2.285, y2 = 3.46 m, V2 = 2.89 m/s, Fr2 = 0.496,
  h_f = 0.625 m, 14%, 61 kW/m, 0.0015 K). Munson et al., 6th ed.:
  Example 2.6 (inclined circular gate, 1.23 MN at 11.6 m, 1.07e5 N m).
  Nothing was fixed. The guide's new "Reference problems" section records
  the one difference of method: Example 6.7's pressure drop is 270 kPa
  with Colebrook's f = 0.0227, against the 265 kPa the book gets from a
  chart-read f and a rounded head loss.

- **B4.materials. Textbook reference problems** (2026-11-02).
  `materials/tests/test_reference_problems.py` checks ten worked examples
  of Callister and Rethwisch, 9th ed., against the printed answers:
  Example 3.3 (copper, 8.89 g/cm^3), 5.1 (carbon flux, 2.4e-9
  kg/(m^2 s)), 5.2 (carburizing, 25,400 s), 5.4 (Cu in Al at 500 C,
  4.2e-14 m^2/s), 6.1 (copper elongation, 0.77 mm), 6.2 (brass rod,
  eps_z = 7.35e-4, F = 5600 N), 8.1 (glass flaw, 8.2 um), Design Example
  8.2 (S-590 rupture, 233 h), 9.2 (Pb-Sn fractions, 0.67 and 0.33) and
  19.1 (brass thermal stress, 106 C). All agree, so nothing was fixed.
  The guide notes the one difference of method: for Example 5.2 the exact
  inverse error function gives 25,370 s (7.05 h) where the book's table
  interpolation prints 25,400 s (7.1 h). Examples from Dowling are left
  for later.

- **B4.heat_transfer. Textbook reference problems** (2026-11-02).
  `heat_transfer/tests/test_reference_problems.py` checks five worked
  examples of Incropera, DeWitt, Bergman and Lavine, 6th ed., against the
  printed answers: Example 1.1 (furnace wall, 1700 W), 1.2 (steam pipe,
  E = 2270 W/m^2, 577 + 421 = 998 W/m), 5.1 (thermocouple junction,
  D = 0.706 mm, Bi = 2.35e-3, 199 C after 5.2 s), 7.4 (Churchill-Bernstein,
  Nu = 40.6 at Re = 6071) and 11.1 (counterflow oil cooler, 8524 W,
  40.2 C, LMTD 43.2 C, L = 65.9 m). All agree, so nothing was fixed or
  added to the guide beyond the test's mention. Example 1.2's irradiation
  is not checked, since 298 K and 298.15 K round it to 447 and 448 W/m^2
  and the printed value could not be confirmed. Examples from Çengel and
  Ghajar are left for later.

- **B4.signal_processing. Textbook reference problems** (2026-11-02).
  `signal_processing/tests/test_reference_problems.py` checks nine worked
  examples against the printed answers. Oppenheim and Schafer, 3rd ed.:
  Example 4.2 (cos 16000πt at 6 kHz aliases to 2 kHz), 7.3 (bilinear
  Butterworth, N = 6, and Ωc = 0.766 meeting the stopband edge), 7.4
  (Butterworth 14, Chebyshev 8, elliptic 6), 7.8 (Kaiser β = 5.653,
  M = 37) and the sec. 7.8.1 Parks-McClellan design (Kaiser's estimate
  M = 26; M = 26 misses δ2 = 0.001 and M = 27 meets it). Proakis and
  Manolakis, 4th ed.: Examples 1.4.3 (aliased frequencies at 5 kHz),
  2.3.2 (linear convolution), 2.6.1 (cross-correlation) and 7.2.1
  (circular convolution). All agree, so nothing was fixed or added to
  the guide; `kaiser_order` uses A − 7.95 where the book uses A − 8,
  which gives the same M here.
- **B4.structural_engineering. Textbook reference problems** (2026-11-02).
  `structural_engineering/tests/test_reference_problems.py` checks six
  worked examples against the printed answers. Hibbeler, *Mechanics of
  Materials*, 10th ed.: Example 8.1 (4 ft cylinder at 20 ksi, 417 psi;
  the sphere 833 psi), 9.2 (element rotated 30° clockwise, -25.8 MPa,
  -4.15 MPa, -68.8 MPa) and 13.1 (pinned W200x46, 4 m, σcr = 320 MPa
  below yield). Gere and Goodno, 7th ed.: Examples 3-1 (solid bar in
  torsion, 4530 psi, 1.62°, T_max = 3980 lb-in), 7-2 (element at 45°)
  and 7-3 (principal stresses 13,500 and -5,440 psi at -14.8°, τmax =
  9,490 psi). All agree, so nothing was fixed or added to the guide.
- **B4.control_systems. Textbook reference problems** (2026-11-02).
  `control_systems/tests/test_reference_problems.py` checks seven worked
  examples of Ogata, *Modern Control Engineering*, 5th ed., against the
  printed answers: Example 5-1 (ζ = 0.6, ωn = 5: tp = 0.785 s,
  Mp = 0.095, ts = 1.33 s), 5-12 (Routh column 1, 2, 1, -6, 5), 6-1 (root
  locus of K/(s(s+1)(s+2)): breakaway -0.4226 at K = 0.3849, ±j√2 at
  K = 6), 7-26 (lead for 4/(s(s+2)), Kv = 20, 50°: α = 0.24, ωc = 9,
  41.7(s+4.41)/(s+18.4)), 8-1 (Ziegler-Nichols, Kcr = 30, Pcr = 2.8099,
  6.3223(s+1.4235)^2/s), 10-1 (K = [199 55 8]) and 10-6 (observer
  Ke = [120.6 20]). Example 10-6 found a bug: `place_poles` rejects a
  pole repeated more often than rank(B), so `state_feedback_poles` and
  `observer_gain` now use Ackermann's formula for one input or output
  with repeated poles. Example 7-26 differs by method: Ogata reads 17°
  and 9 rad/s off the Bode plot (17.96° and 8.95 rad/s exactly), which
  moves his zero and pole about 0.05 rad/s; the guide records it. The
  example numbers were not checked against a printed copy (5-12 and 10-6
  least sure). Examples from Franklin, Powell and Emami-Naeini are left
  for later.

- **B4.thermodynamics. Textbook reference problems** (2026-11-02).
  `thermodynamics/tests/test_reference_problems.py` checks ten worked
  examples against the printed answers. Çengel and Boles, 9th ed.:
  Examples 10-1 (26.0%, x = 0.886), 10-3a, b and c (33.4%, 37.3%, 43.0%;
  x = 0.813, 0.915, 0.804), 10-4 (reheat, 45.0%, x = 0.896), 10-5 (open
  heater, 46.3%, y = 0.227), 9-5 (cold-air Brayton, 0.448) and 15-8
  (octane with 400% air, 962 K). Moran et al., 9th ed.: Examples 8.1
  (37.1%, x = 0.675, 3.77e5 kg/h) and 8.2 (31.4%). Nothing was fixed.
  Two quantities meet one digit fewer than printed because the book's
  steam tables differ from IF97 (y = 0.2271 against 0.2270, Moran's
  x = 0.6748 against 0.6745). The guide's new "Textbook reference
  problems" section records two differences of method, each bounded by
  a test: Example 11-1's Peng-Robinson R-134a gives COP 3.94 against
  3.97, and Example 9-6's constant cp gives a 805 K turbine exit against
  the air tables' 853 K (r_bw 0.615 against 0.592, 27.4% against 26.6%).
  The Chapter 10 answers were confirmed against a published
  reproduction of the examples; the rest were not checked against a
  printed copy (9-6 and 15-8 least sure).

- **B4.mechanical_engineering. Textbook reference problems** (2026-11-02).
  `mechanical_engineering/tests/test_reference_problems.py` checks eight
  worked problems. Budynas and Nisbett, *Shigley's Mechanical Engineering
  Design*, 10th ed.: Example 14-1 (Lewis rating, 365 lbf, 6.95 hp) and
  the solutions-manual answers to Problems 7-1 (DE-Goodman shaft,
  27.27 mm), 10-3 and 10-5 (spring rate, Bergsträsser factor, load and
  stress at solid) and 11-8 to 11-10 (bearing catalog rating). Rao,
  *Mechanical Vibrations*, 5th ed.: Example 6.11 (three-mass chain).
  Nothing was fixed. The guide's new "Textbook reference problems"
  section records Shigley's rounded intermediates (K_v = 1.52 gives
  365 lbf, unrounded 364 lbf) and Rao's low-precision roots of the
  cubic (1.2471 and 1.8025 against 1.24698 and 1.80194). The chapter 7,
  10 and 11 answers were read from the published solutions manual; the
  Example 14-1 answers from a reproduction of the example, and Example
  6.11 from a copy of Rao's chapter slides, not from a printed book.
  Shigley's in-text Examples 8-1 (power screw), 8-2 (frustum stiffness)
  and 10-1 could not be checked against a copy and are left for later.

- **B4.electrical_engineering. Textbook reference problems** (2026-11-02).
  `electrical_engineering/tests/test_reference_problems.py` checks six
  worked examples from Glover, Overbye and Sarma, *Power System Analysis
  and Design*: 6th ed. Examples 9.3 (bolted single-line-to-ground fault,
  I1 = -j1.964 pu) and 11.1 (swing angle 16.75 degrees three cycles after
  a terminal fault), and 7th ed. (with Birchfield, 2022) Examples 6.9
  (five-bus Newton-Raphson load flow: Y22, all bus voltages, slack and
  PV-bus generation), 12.4 (maximum swing 0.7003 rad) and 12.5 (critical
  clearing time just under 0.19 s). Nothing was fixed. Example 12.4's
  0.7003 is 1e-4 above the root of the book's own equation (0.70020;
  `swing_response` gives 0.70025), and the guide's new "Textbook
  reference problems" section says so. The answers were read from
  Overbye's ECEN 460 slides (2019, 2025), Lee's Baylor ELC 4340 slides
  (2022) and two reproductions of Example 11.1, not from a printed book.
  Hayt, Kemmerly and Durbin, the circuits text named in D4, could not be
  checked against a copy, so the DC, transient, op-amp, diode and machine
  models have no reference problems yet.

- **B5. ValidityWarning, piloted in heat_transfer and fluid_engineering**
  (2026-11-01). `_validation.py` has `ValidityWarning(UserWarning)`,
  exported as `engineeringkit.ValidityWarning`, and `validity(model,
  name, value, low, high)`, which warns with the stated range and the
  first offending value and leaves the result unchanged. Eleven checks in
  ten functions: `dittus_boelter_nusselt` (Re >= 1e4, 0.6 <= Pr <= 160),
  `flat_plate_nusselt` (0.6 <= Pr <= 60), `churchill_bernstein_nusselt`
  (Re Pr >= 0.2), `zukauskas_nusselt` (0.7 <= Pr <= 500; ST/SL >= 0.7 for
  aligned banks), `hausen_nusselt` (Re <= 2300), laminar
  `churchill_chu_nusselt` (Ra <= 1e9), `nusselt_film_condensation` (film
  Reynolds number <= 1800), `stokes_settling_velocity` (Re of the result
  <= 1), `rehbock_discharge_coefficient` (H/P <= 5) and `pipe_wave_speed`
  (D/e >= 25). Colebrook in laminar flow, named in the item, needed
  nothing: it returns 64/Re by design, and the transition range already
  raises. `lumped_temperature` and `chilton_colburn_mass_transfer` stay
  unchecked, since their limits need a conductivity or a Prandtl number
  they are not given. The new warnings found four gallery examples
  outside a stated range: the Pohlhausen comparison now stops at Pr = 60,
  the hot-wire anemometer starts at 1.5 m/s (Re Pr = 0.2), the PVC pipe
  in the water-hammer example has a 12 mm wall (D/e = 25), and the
  Stokes example silences the warning where the mismatch is its point.
  The other ten domains have B5.<domain> lines under B5.

- **B2. Dependency floors and Windows tested** (2026-10-31). CI has a
  `test-min-deps` job on Python 3.10 with NumPy 1.24, SciPy 1.10 and
  Matplotlib 3.6 pinned (and pyparsing 3.0, which matplotlib 3.6 needs
  under warnings-as-errors), and `windows-latest` joins the test matrix.
  The floors hold: only two tests and two gallery examples called
  `np.trapezoid` (NumPy 2.0) and now use SciPy's; all tests and all 327
  examples pass locally at the floors. For Windows, four test
  `read_text()` calls on UTF-8 pages now pass `encoding="utf-8"` (cp1252
  is Windows' default), and `MPLBACKEND` moved to the job's `env`, since
  PowerShell rejects `VAR=x cmd`. The Windows row has not run yet: it
  first runs on the next push.

- **B1. Branch coverage** (2026-10-31). `branch = true` is set, so the
  CI coverage job's `fail_under = 100` now fails on a partial branch. By
  0.2.0 the 44 partials had grown to 72, all covered with real tests and
  no pragmas. Most were plots never handed an existing axes or drawn
  without a label. Four found or forced code changes: the Neuber notch
  solver stalled at its 100-step cap for large n (n = 200 returned
  σε = 6e82 instead of 8e5) and now runs Newton on ln σ until it
  converges; `method_of_joints` names improper (parallel or concurrent)
  supports instead of reporting a stall; and the Newton-Raphson load
  flow and moment distribution loops no longer have an exit only `break`
  could take.

- **A2. No generated plot images in the source tree** (2026-10-31).
  The current build (matplotlib 3.11 plot directive) writes `.. plot::`
  figures only under `docs/build/plot_directive/` and `_images/`; neither
  `make html` nor `make html-noplot` puts anything in `docs/source/`. The
  14 `<page>-N.png` files committed in c953e721 were stale leftovers from
  an earlier build that no page referenced, so they are deleted, and
  `.gitignore` ignores `<page>-N.png` in `exercises/` and `tutorials/` in
  case an older toolchain writes them again.

- **A4. Stale release and verification text removed** (2026-10-31).
  `VERIFICATION.md` (0.1.0-era numbers) is deleted, and its check list,
  with today's numbers, is "Reproducing the checks" in
  `docs/source/validation.rst`. `limitations.rst` no longer says nothing
  is on PyPI, and the README's "Release scope (v0.1)" is now "Scope and
  checks".

- **A5. Unreleased section in CHANGELOG.md** (2026-10-31). The
  `## [Unreleased]` heading went in with 0.2.0. The two commits the item
  named (60ee4ee5, dated gallery example titles; 1f8cda4b, the Linux-only
  control exercise check) shipped in 0.2.0, so their lines are under
  `## [0.2.0]`, not Unreleased.

### 0.2.0 (2026-10-30)

- **A1. Pending working tree committed** (2026-10-30). `docs/Makefile`
  and `docs/make.bat` went in with the `docs.yml` and `release.yml`
  changes that call `cd docs && make html` (428268e9), and the 14
  exercise and tutorial plot PNGs followed in c953e721. `git status` is
  clean and the Docs workflow is green on `main` (run 41043990280).

- **G1. Release 0.2.0** (2026-10-30). Tagged `v0.2.0`; PyPI, the GitHub
  Release and the docs rebuild came from `release.yml`, and Zenodo minted
  the version DOI 10.5281/zenodo.23795071, now in `CITATION.cff` and the
  README badge.

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
