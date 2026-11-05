# JerkSimulator

![Version: V4.3.2](https://img.shields.io/badge/version-V4.3.2-blue?style=plastic&labelColor=555) ![HTML5: single-file](https://img.shields.io/badge/HTML5-single--file-orange?style=plastic&labelColor=555) ![Offline: ready](https://img.shields.io/badge/offline-ready-green?style=plastic&labelColor=555) [![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue?style=plastic&labelColor=555)](LICENSE)

JerkSimulator is a single-file HTML5 simulator for validated motion profiles. It supports jerk-limited S-curves, classic linear ramps, and hybrid profiles for industrial motion analysis.

## Live application:
[https://firerayo.github.io/JerkSimulator/](https://firerayo.github.io/JerkSimulator/)

## Features

- Calculate exact seven-segment jerk-limited S-curve profiles without Euler integration.
- Use **Single jerk** mode or independent `J1`, `J2`, `J3`, and `J4` limits.
- Set a jerk limit to `0` to disable S-curve smoothing for that transition.
- Model classic linear velocity ramps with instantaneous acceleration changes.
- Show ideal `+∞` and `-∞` jerk impulses for instantaneous acceleration steps.
- Combine S-curve and linear-ramp transitions in hybrid profiles.
- Configure target distance, maximum speed, acceleration, and deceleration.
- Configure optional initial speed `V0` and final speed `Vf`.
- Automatically select `Vmax + cruise`, `Vmax without cruise`, or a reduced peak speed.
- Validate distance, speed, acceleration, deceleration, jerk, and final solver residuals.
- Show synchronized Position, Velocity, Acceleration, and Jerk graphs.
- Inspect graphs with a synchronized cursor, zoom, pan, and segment boundaries.
- Overlay the previous simulation for direct comparison.
- Animate the exact calculated position `x(t)`.
- Adjust animation speed and repeat the animation.
- Export high-precision CSV data with profile metadata.
- Print the interface and graphs or save them as PDF through the browser print dialog.
- Use English, Spanish, Italian, or Brazilian Portuguese.
- Select Light or Dark Neomorphism themes.
- View the complete GNU GPLv3 license inside the application without Internet access.

## Quick start

1. Download `index.html`.
2. Open it in a modern browser.
3. Enter the required motion limits.
4. Select **Single jerk** or **Independent jerks**.
5. Select **Simulate**.
6. Inspect the calculated values and graphs.
7. Use **Export CSV** or **Print to PDF** when required.

The application starts with this default profile:

```text
Target distance S = 2000 mm
Max speed Vmax = 500 mm/s
Max acceleration Amax = 300 mm/s²
Max deceleration Dmax = 300 mm/s²
Max jerk j = 500 mm/s³
Initial speed V0 = 0 mm/s
Final speed Vf = 0 mm/s
Plot / CSV samples = 1000
```

The default simulation is generated automatically when the application opens. **Boundary velocities and output sampling** remains collapsed until you open it.

## Compatibility and limitations

JerkSimulator uses HTML5, CSS3, Vanilla JavaScript, and Canvas. It has no required runtime framework, package manager, build process, CDN, or web service.

The application is designed to work directly from `file://` in modern browsers. Browser printing and PDF output can vary by browser and operating system.

The motion model is ideal kinematics. It does not model motor torque, load inertia, backlash, elasticity, friction, resonance, or control-loop dynamics. Always validate final parameters on the real machine and against the drive documentation.

## Data and privacy

Motion calculations run locally in the browser. The main application functions do not require a network connection.

The application stores only language and theme preferences in `localStorage` when that storage is available. Simulation inputs and generated CSV data remain local to the browser.

GitHub links require a network connection only when you open them.

## Development and verification

The readable `index.html` file is the application source. No compilation step is required.

For V4.3.2, static verification confirmed no external runtime script or stylesheet dependencies. The executable JavaScript passed a syntax check with `node --check`.

The public GitHub Pages site and repository can contain an earlier deployed version until the updated `index.html` is published.

## Contributing and support

Repository:
https://github.com/FireRayo/JerkSimulator/

Issues:
https://github.com/FireRayo/JerkSimulator/issues

Use GitHub Issues to report reproducible problems or propose changes.

## License and credits

Copyright (C) 2026 Raymundo Ortiz.

JerkSimulator is licensed under the **GNU General Public License, version 3 or, at your option, any later version** (`GPL-3.0-or-later`).

See [LICENSE](LICENSE) for the complete GNU GPLv3 text.

This program is distributed without warranty, as specified by the GNU General Public License.

The V4.3.2 source identifies no included third-party components. Third-party components added in future versions must retain their original copyright notices and licenses.
