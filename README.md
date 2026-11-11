# Assignment 5: Vibrating Box

Run `test_box` to check the spring, rigid body transform, and RK simulation. Run
`simulate_box` to integrate the example system and write `vibrating_box.mp4` in
the repository root.

## Vibration modes

In MATLAB, from this directory:

```matlab
test_box_modes;           % Existing tests, solvers, modal checks, short video test
results = analyze_box_modes();  % Three report figures and the selected 1080p video
verify_mode_artifacts;   % Check saved figures, video metadata, and decoded frames
```

To regenerate the analysis and figures without rendering the full video:

```matlab
results = analyze_box_modes(struct('make_video',false));
```

The modal driver writes to the repository root by default. It does not
overwrite the original `vibrating_box.mp4`. All MATLAB source and tests also
live in the root. Tests use a temporary directory. No additional MATLAB
toolboxes are required.

Repository deliverables, all in the root:

- All MATLAB source/tests (`.m`) and this `README.md`.
- `mode_1_comparison.png` through `mode_3_comparison.png`: 300 dpi, stacked
  x, y, and orientation plots comparing the analytical cosine and nonlinear RK4.

Local artifacts (kept on disk, not tracked):

- Matching `.fig` files: editable MATLAB figures.
- `mode_1_nonlinear.mp4`: selected actual nonlinear trajectory, 1920 x 1080,
  60 fps, approximately 24 seconds. Mode selection uses the largest validated
  corner travel. The current parameters select mode 1.
- `vibrating_box.mp4`: original example animation.
- `validation.md`: numerical results and metric definitions.
- `modal_results.mat`: cached parameters, equilibrium, Jacobian, mass-normalized
  mode shapes, frequencies, amplitude trials, trajectories, and validation metrics.
  `verify_mode_artifacts` reads this cache from the root.
- `validation.json`: machine-readable numerical results.
- `artifact_validation.json` and `decoded_frame_*.png`: saved-media checks and
  sampled frames decoded from the MP4.

Extra artifacts are excluded locally through `.git/info/exclude`, including
all videos, editable `.fig` files, `validation.md`, logs, `.mat` caches, JSON
diagnostics, and sampled frames. The remaining `outputs/` directory is retained
locally and excluded as well. No committed `.gitignore` is used. The existing
`Ass5.pdf` exclusion is preserved, and `.DS_Store` stays local but is not tracked.
Nothing is deleted; only MATLAB files, the three comparison PNGs, and README
are intended for Git.

The analysis cache and JSON diagnostics can be regenerated with
`analyze_box_modes`. The figure-only command above does not render a video or
rebuild its saved metadata; to verify an existing saved video, retain its
matching cache or run the full analysis first.

The equilibrium comes from the Assignment 3 Newton solver. The rate-function
Jacobian uses the Assignment 3 central-difference routine. Modes come from
`Q = -A(4:6,1:3)` and `eig(Q)`, sorted by increasing frequency and normalized
with `Phi'*M*Phi = I`.

Each trajectory starts with `Veq + [epsilon*Phi(:,mode); zeros(3,1)]`.
The amplitude search chooses the largest tested corner displacement with
mass-weighted RMS error, peak error, and unwanted-mode leakage all below 1%
over 12 predicted periods. Per-coordinate errors are reported separately and
can exceed 1% for small mode-shape components. The nonlinear trajectory is
predominantly one mode, not an exactly single-frequency solution.

Validation includes equilibrium residual and rest simulation, Jacobian block
structure and independent five-point reference, mass-weighted stiffness
symmetry, positive eigenvalues, eigenpair residuals, energy drift, measured
periods, and fourth-order RK4 convergence under step halving. The finest RK4
trajectory is used in every exported artifact.

The animation uses fixed equal-scale axes, an external legend, and a fixed
attachment close-up in true coordinates. Playback is slowed to two seconds
per predicted period, with the slowdown labeled. There is no motion gain.
Time and position follow the assignment's dimensionless convention `(-)`;
orientation is in `(rad)`, angular frequency in `rad/(-)`, and frequency in
`cycles/(-)`, not Hz.

`animate_box(t,V,params)` keeps the original output name and 30 fps playback.
An optional fourth struct configures `output_file`, `fps`, `figure_size`,
`playback_speed`, `mode`, `omega`, `equilibrium`, `legend_location`, `visible`,
`sample_frames`, and `focus_inset`. See `analyze_box_modes` for an example.

YouTube uploading is a separate authenticated browser step, not part of the
MATLAB scripts. Only the selected MP4 is intended for an Unlisted upload.
Completed Unlisted upload: https://youtu.be/GTtZjoMQH-o.
Upload verification is retained locally in the ignored
`outputs/modes/youtube_upload.md`.
The unavailable video example and Orion reference were not reviewed.

### Source attribution

`approximate_jacobian.m` and `multi_newton_solver.m` are copied unchanged except
for attribution headers from this user's `../MTH3199-Assignment-3`, commit
`2c614d9f8d148b162ed373f79ef71a0706bce013`. That repository is not modified.
The shared parameters retain the original values; `l0_list` is now explicitly
a row vector rather than an accidentally square constant array.

## Files

- `compute_spring_force.m`: spring force from its endpoints, stiffness, and natural length.
- `compute_rbt.m`: transforms box-frame points into world coordinates.
- `compute_accel.m`: sums spring forces and moments, then applies linear and angular momentum.
- `box_rate_func.m`: first-order state derivative for `[x;y;theta;vx;vy;vtheta]`.
- `explicit_RK_fixed_step_integration.m`: fixed-step explicit Runge-Kutta integrator and RK step.
- `simulate_box.m`: example parameters, initial state, integration, and animation call.
- `animate_box.m`: fixed-axis spring and box animation saved as MPEG-4.
- `test_box.m`: basic numerical checks.
- `box_parameters.m` and `rk4_tableau.m`: shared geometry and integration method.
- `analyze_box_modes.m`: equilibrium, modes, amplitude search, and convergence.
- `box_energy.m`: total mechanical energy for numerical checks.
- `export_mode_comparison.m`: LaTeX report figures and editable `.fig` exports.
- `write_modal_validation.m`: machine-readable results and numerical report.
- `test_box_modes.m` and `verify_mode_artifacts.m`: numerical and media checks.
