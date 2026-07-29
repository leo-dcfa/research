# State vectors in space flight

The **state vector** is the minimal set of numbers that, together with a dynamics model and an epoch, determines a spacecraft's entire past and future trajectory. In astrodynamics it is almost always the **6-vector** of inertial position and velocity at a time $t$:

$$\mathbf{x}(t) = \begin{bmatrix} \mathbf{r} \\ \mathbf{v} \end{bmatrix} = \begin{bmatrix} x & y & z & \dot x & \dot y & \dot z \end{bmatrix}^{\mathsf T}$$

- **Why 6?** Motion obeys a second-order ODE $\ddot{\mathbf{r}} = \mathbf{a}(\mathbf{r}, \dot{\mathbf{r}}, t)$ in 3D. Rewriting it as a first-order system needs 3 positions + 3 velocities → the state space is $\mathbb{R}^6$. Six initial conditions, six constants of motion.
- A state vector is **meaningless without three tags**: the **epoch** (when), the **frame** (relative to what axes/origin), and the **units** (usually km and km/s).
- It's the *lingua franca* of flight dynamics: JPL Horizons, SPICE kernels, GPS broadcast ephemerides, and every onboard navigation filter speak in state vectors.

## The equation of motion

Under a central body plus perturbations, the state obeys

$$\dot{\mathbf{x}} = f(\mathbf{x}, t) = \begin{bmatrix} \mathbf{v} \\[4pt] -\dfrac{\mu}{r^{3}}\,\mathbf{r} + \mathbf{a}_{p} \end{bmatrix}, \qquad r = \lVert \mathbf{r} \rVert$$

with $\mu = GM$ the gravitational parameter ($\mu_\oplus = 398\,600.4418\ \mathrm{km^3/s^2}$) and $\mathbf{a}_p$ everything else. Set $\mathbf{a}_p = 0$ and you have the **two-body problem**, which is integrable in closed form.

Constants of the unperturbed motion, all computable directly from the state:

| Quantity | Formula | Meaning |
|---|---|---|
| Specific angular momentum | $\mathbf{h} = \mathbf{r} \times \mathbf{v}$ | Fixes the orbit plane |
| Specific energy | $\xi = \dfrac{v^{2}}{2} - \dfrac{\mu}{r}$ | Sets the size: $a = -\mu / (2\xi)$ |
| Eccentricity vector | $\mathbf{e} = \dfrac{\mathbf{v} \times \mathbf{h}}{\mu} - \dfrac{\mathbf{r}}{r}$ | Points at periapsis, magnitude $e$ |
| Vis-viva | $v^{2} = \mu\left(\dfrac{2}{r} - \dfrac{1}{a}\right)$ | Speed anywhere on the orbit |

- $\xi < 0$ → ellipse (bound), $\xi = 0$ → parabola (escape), $\xi > 0$ → hyperbola (flyby / interplanetary).
- Bound orbits also give the period $T = 2\pi\sqrt{a^{3}/\mu}$ and mean motion $n = \sqrt{\mu/a^{3}}$.

## Frames — the part that actually bites

The numbers are only as good as the axes they're expressed in. Confusing frames is a classic mission-losing bug.

| Frame | Origin / axes | Used for |
|---|---|---|
| **ECI / GCRF (J2000, ICRF)** | Earth centre, inertial, non-rotating w.r.t. distant quasars | Orbit propagation, integration — Newton's laws only hold here |
| **ECEF / ITRF** | Earth centre, rotates with the Earth | Ground tracks, geodetic lat/lon/alt, GNSS products |
| **Heliocentric (HCI / EMO2000)** | Sun centre | Interplanetary cruise |
| **RIC / RSW (LVLH)** | Spacecraft-centred: radial, in-track, cross-track | Covariance reporting, conjunction assessment, formation flying |
| **Perifocal (PQW)** | Orbit plane: $\hat{\mathbf{p}}$ to periapsis | Intermediate step in element ↔ state conversions |

- ECI → ECEF is not a single rotation: it's **precession → nutation → Earth rotation angle → polar motion** (IAU-76/FK5 or IAU-2000/2006 CIO chains).
- **Time scales matter too**: dynamics integrate in **TDB/TT**, Earth rotation needs **UT1**, telemetry is stamped in **UTC** (leap seconds!), GNSS uses **GPS time** (no leap seconds, currently $\mathrm{TAI} - 19\ \mathrm{s}$). A 1 s error in UT1 smears the ground track by ~465 m at the equator.
- **Inertial ≠ non-accelerating**: ECI orbits the Sun, but for near-Earth work the tidal residual is negligible, so it is treated as inertial.

## The other chart: Keplerian elements

The same six numbers can be written as **classical orbital elements** (COE) — a different coordinate chart on the same 6D state space, chosen so that five of the six are *constant* in two-body motion.

| Element | Symbol | Role |
|---|---|---|
| Semi-major axis | $a$ | size |
| Eccentricity | $e$ | shape |
| Inclination | $i$ | tilt of orbit plane |
| RAAN | $\Omega$ | where the plane cuts the equator |
| Argument of periapsis | $\omega$ | orientation of the ellipse in-plane |
| True anomaly | $\nu$ | **where the spacecraft is right now** (the only fast variable) |

- **State ↔ elements is an exact, invertible map** (see `state_vectors.py`), except at degeneracies: circular orbits ($e \to 0$) lose $\omega$, equatorial orbits ($i \to 0$) lose $\Omega$. Fixes: use equinoctial elements, or the singularity-free workarounds (argument of latitude $u = \omega + \nu$, true longitude).
- Anomaly zoo: **true** $\nu$ (geometric angle) $\to$ **eccentric** $E$ $\to$ **mean** $M = E - e\sin E$ (linear in time, $M = M_0 + n\,\Delta t$). Going from time to position means inverting **Kepler's equation** — transcendental, solved by Newton iteration.
- **State vectors are better for**: numerical integration, maneuvers, filters, anything with perturbations. **Elements are better for**: mission design intuition, secular drift analysis, catalogue storage.

## Propagation — moving the state through time

**1. Analytic (two-body only).** Solve Kepler's equation, then rebuild the state with the **Lagrange $f$ and $g$ coefficients**, which express the future state as a linear combination of the current one *in the orbit plane*:

$$\mathbf{r}(t) = f\,\mathbf{r}_0 + g\,\mathbf{v}_0, \qquad \mathbf{v}(t) = \dot f\,\mathbf{r}_0 + \dot g\,\mathbf{v}_0, \qquad f\dot g - \dot f g = 1$$

Universal-variable formulations make this work for ellipse, parabola and hyperbola with one code path.

**2. Numerical (real missions).** Integrate $\dot{\mathbf{x}} = f(\mathbf{x},t)$ with a high-order integrator (RKF7(8), Dormand–Prince, Gauss–Jackson for long arcs). Perturbations, in rough order of magnitude in LEO:

| Perturbation | Typical magnitude (LEO) | Notes |
|---|---|---|
| Two-body $\mu/r^{2}$ | $\sim 8.7\ \mathrm{m/s^2}$ | reference |
| $J_2$ oblateness | $\sim 10^{-2}\ \mathrm{m/s^2}$ | dominant perturbation; secular, predictable |
| Atmospheric drag | $10^{-6}$–$10^{-4}\ \mathrm{m/s^2}$ | worst-modelled term; density uncertainty is huge |
| Higher geopotential ($J_3$, tesserals) | $\sim 10^{-5}\ \mathrm{m/s^2}$ | EGM2008, GGM03 models |
| Third body (Moon, Sun) | $\sim 10^{-6}\ \mathrm{m/s^2}$ | dominates in GEO and beyond |
| Solar radiation pressure | $\sim 10^{-7}\ \mathrm{m/s^2}$ | needs area, reflectivity, eclipse model |

Formulations: **Cowell** (integrate total acceleration — simplest, what everyone does now), **Encke** (integrate only the deviation from a reference conic), **variation of parameters** (integrate the *elements* instead, since they change slowly).

**3. The $J_2$ secular rates** — the single most useful closed-form perturbation result:

$$\dot\Omega = -\tfrac{3}{2}\, J_2 \left(\frac{R_\oplus}{p}\right)^{2} n \cos i, \qquad \dot\omega = \tfrac{3}{4}\, J_2 \left(\frac{R_\oplus}{p}\right)^{2} n \left(5\cos^{2} i - 1\right)$$

with $p = a(1-e^2)$ and $J_2 \approx 1.0826 \times 10^{-3}$.

- Choose $i \approx 98^\circ$ so $\dot\Omega \approx 0.9856^\circ/\text{day}$ and the plane tracks the Sun → **sun-synchronous orbit**.
- Choose $i = 63.4^\circ$ so $5\cos^2 i - 1 = 0$ and periapsis stops drifting → **Molniya orbit**.

**4. Osculating vs mean.** The instantaneous state converted to elements gives **osculating** elements (wobbling every orbit). Averaging out the short-period terms gives **mean** elements — smoother, and what **TLEs** contain. A TLE is *not* a state vector: it is a set of Brouwer mean elements that is only valid when fed to **SGP4**, which bakes in its own perturbation model.

## Uncertainty: the state vector is never exact

Real navigation carries a **state estimate plus a covariance**:

$$\hat{\mathbf{x}} \in \mathbb{R}^{6}, \qquad P = \mathbb{E}\!\left[(\mathbf{x}-\hat{\mathbf{x}})(\mathbf{x}-\hat{\mathbf{x}})^{\mathsf T}\right] \in \mathbb{R}^{6\times 6}$$

- Uncertainty is moved forward by the **state transition matrix** $\Phi(t, t_0) = \partial \mathbf{x}(t) / \partial \mathbf{x}(t_0)$, which itself satisfies the variational equation

$$\dot\Phi = A(t)\,\Phi, \qquad A = \frac{\partial f}{\partial \mathbf{x}} = \begin{bmatrix} 0 & I \\[4pt] \dfrac{\partial \mathbf{a}}{\partial \mathbf{r}} & \dfrac{\partial \mathbf{a}}{\partial \mathbf{v}} \end{bmatrix}, \qquad \Phi(t_0,t_0) = I$$

  For two-body gravity the gravity gradient is $\dfrac{\partial \mathbf{a}}{\partial \mathbf{r}} = -\dfrac{\mu}{r^{3}}\left(I - 3\,\hat{\mathbf{r}}\,\hat{\mathbf{r}}^{\mathsf T}\right)$ — the tidal tensor.
- Covariance then propagates as $P^{-} = \Phi P \Phi^{\mathsf T} + Q$, with $Q$ the process noise absorbing unmodelled drag, SRP and thrust.
- **Error ellipsoids stretch along-track.** A small semi-major-axis error becomes a period error, which integrates into a growing in-track offset ($\propto \Delta t$). This is why covariance is reported in the RIC frame and why conjunction screening cares mostly about the in-track term.
- **Orbit determination** fits the 6 (or more) state parameters to observations — range, Doppler/range-rate, angles, GNSS pseudoranges — by weighted least squares (batch) or a Kalman/EKF/UKF (sequential). Typical filter state is augmented well past 6: drag coefficient, SRP coefficient, clock bias and drift, empirical accelerations, station biases.

## Maneuvers and relative motion

- An **impulsive burn** is the cleanest thing in the whole model: it leaves $\mathbf{r}$ untouched and adds to the velocity block,

$$\mathbf{x}^{+} = \mathbf{x}^{-} + \begin{bmatrix} \mathbf{0} \\ \Delta\mathbf{v} \end{bmatrix}$$

  Everything about the new orbit follows by re-deriving the elements from the new state. Finite burns are modelled by adding a thrust term to $\mathbf{a}_p$ and integrating through.
- **Lambert's problem** — the transfer workhorse: given $\mathbf{r}_1$, $\mathbf{r}_2$ and a time of flight, find the connecting conic. Differencing its endpoint velocities against the departure and arrival states gives the two $\Delta v$'s; sweeping departure and arrival dates gives the **porkchop plot** used for every interplanetary launch window.
- **Rendezvous** uses the *relative* state $\delta\mathbf{x}$ in the target's LVLH frame. Linearised about a circular orbit it gives the **Clohessy–Wiltshire** equations,

$$\ddot x - 3n^{2}x - 2n\dot y = 0, \qquad \ddot y + 2n\dot x = 0, \qquad \ddot z + n^{2}z = 0$$

  whose closed-form STM is what drives approach-corridor planning (and explains the counter-intuitive result that thrusting forward makes you fall *behind*).

## Practical notes

- **SPICE** (NAIF) is the standard toolkit: SPK kernels store ephemerides, `spkezr` returns a 6-element state in a requested frame with light-time correction. Interplanetary work is essentially all SPICE-mediated state vectors.
- Standard interchange formats: **OEM/OMM** (CCSDS) for state vectors and mean elements, **SP3** for GNSS precise orbits, **TLE/3LE** for the public catalogue.
- Watch for: km vs m, degrees vs radians, position–velocity ordering, frame epoch (J2000 vs of-date), and whether a "state" already has light-time/aberration corrections applied.
- Sanity checks that catch most bugs: recompute $\xi$ and $\lVert \mathbf{h}\rVert$ before and after a propagation (they must be conserved in two-body), and round-trip state → elements → state.

## Code

- `state_vectors.py` — state ↔ classical elements both ways, invariants, an RK4 propagator with $J_2$, and a check of the predicted nodal regression against the integrated orbit.
- `covariance_stm.py` — integrate the variational equation for $\Phi$, validate it against finite differences, and propagate a covariance to watch the in-track error blow up.

## References

- Vallado, *Fundamentals of Astrodynamics and Applications* — the standard reference, algorithms given in pseudocode.
- Curtis, *Orbital Mechanics for Engineering Students* — gentler derivations of the conversions above.
- Tapley, Schutz & Born, *Statistical Orbit Determination* — the covariance/STM/filtering half.
- Montenbruck & Gill, *Satellite Orbits* — models and numerical methods.
