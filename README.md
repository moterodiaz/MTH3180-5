# Assignment 5: Vibrating Box

Run `test_box` to check the spring, rigid body transform, and RK simulation. Run
`simulate_box` to integrate the example system and write `vibrating_box.mp4` in
the repository root.

## Files

- `compute_spring_force.m`: spring force from its endpoints, stiffness, and natural length.
- `compute_rbt.m`: transforms box-frame points into world coordinates.
- `compute_accel.m`: sums spring forces and moments, then applies linear and angular momentum.
- `box_rate_func.m`: first-order state derivative for `[x;y;theta;vx;vy;vtheta]`.
- `explicit_RK_fixed_step_integration.m`: fixed-step explicit Runge-Kutta integrator and RK step.
- `simulate_box.m`: example parameters, initial state, integration, and animation call.
- `animate_box.m`: fixed-axis spring and box animation saved as MPEG-4.
- `test_box.m`: basic numerical checks.
