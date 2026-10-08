# Model notes

These notes are here so a reader can follow the model without guessing what the code means. The project is a numerical experiment, not an operational weather model.

## The experiment in one sentence

Run the same small atmospheric toy model in `float` and `double`, then change one starting value by one representable step at a time in each precision and measure how much the outcome changes.

## What the variables mean

Each of the 20 sites has four values:

- `circulation` is the Lorenz-96 model variable. It represents a simplified atmospheric state; it is not wind speed or temperature in real-world units.
- `humidity` is a dimensionless moisture proxy.
- `cloud` is a dimensionless cloud-water proxy.
- `rain` is accumulated rain proxy for the whole run.

The sites form a ring. For example, the site after site 19 is site 0. This boundary choice is part of the Lorenz-96 test model, not a map of a real region.

## The atmospheric motion equation

For site `i`, Lorenz-96 uses:

```text
dx[i]/dt = (x[i+1] - x[i-2]) * x[i-1] - x[i] + F
```

The indices wrap around the ring, and `F` is the steady forcing. This compact equation produces chaotic motion that is useful for studying how small changes spread. It leaves out the full three-dimensional fluid equations used in serious weather prediction.

## The moisture and rain experiment

The moisture, cloud, and rain terms are intentionally simple additions for this project. They are not claimed to represent real cloud microphysics.

Condensation is allowed when the moisture proxy is above a saturation threshold and circulation is above a lift threshold:

```text
condensation = rate * max(0, humidity - saturation) * max(0, circulation - lift threshold)
```

Humidity relaxes toward a fixed background value, diffuses between neighboring sites, and decreases by condensation. Cloud proxy increases through condensation and decreases through dissipation and rain. Rain accumulates when cloud proxy passes a threshold. The constants are collected at the top of `src/main.cpp` so a reader can find and change them in one place.

## The one-ULP perturbation

`std::nextafter(value, +infinity)` returns the next representable value of the chosen type above `value`. For a value near 8, one `float` step is larger than one `double` step. Stormline makes this change to the initial circulation at one site, then reruns the same simulation with the same time step. The CSV also compares the unperturbed `float` and `double` baselines.

This is a controlled perturbation experiment. It asks whether a tiny initial difference can grow in this model and how the result changes between two arithmetic types. It does not claim that a particular computer's hardware made an error. The one-ULP probe is applied only at the start of each run; the next research step would study precision effects during the calculations themselves.

## Time integration

The program uses classical fourth-order Runge–Kutta (RK4). At each time step, it evaluates the rate of change four times—at the beginning, twice in the middle, and at the end—then combines those rates with weights `1, 2, 2, 1`. This is more accurate for a smooth trajectory than taking one forward-Euler step of the same size, at the cost of four rate evaluations per step.

## What the CSV measures

The first CSV row compares the unperturbed `float` and `double` baselines. The remaining rows compare each one-ULP run with its baseline in the same precision:

- `ulp_delta` records the actual size of the representable step at the nudged value; it is different for `float` and `double`.
- `baseline_storm_count` and `compared_storm_count` count sites above the rain threshold in each run.
- `changed_storm_site_count` and `changed_storm_site_ids` show the sites whose threshold status changed; IDs are separated by `|`.
- `steps` and `time_step` record the settings used for that row.
- `baseline_rain_total` and `compared_rain_total` give the total rain proxy in each run.
- `rainfall_absolute_difference` adds up the absolute rain-proxy differences over all sites.
- `final_state_rms_difference` measures the end-of-run separation in circulation, humidity, and cloud values.

These metrics describe this model and these settings. They are not real-world forecast skill scores.

## The sensitivity atlas

`stormline.svg` shows a local final-state difference for every affected-site/nudged-site pair. The left panel is `double`; the right panel is `float`. Both panels use the same color scale, so their shading can be compared directly. An orange cell border means the rain-threshold status changed at that affected site for that particular nudge. The atlas helps answer two different questions: how far a small change spread, and whether it changed the simple storm-like classification.

## What would make the next version stronger

1. Add output for the baseline rain and storm-like state at every site, so a reader can inspect where the threshold was crossed.
2. Repeat the experiment over several time steps and model sizes to see whether the result depends on the numerical setup.
3. Probe rounding during the calculation, not only in the initial condition, and see whether the most influential sites change.
4. Only after profiling, consider an assembly implementation of a measured hot loop. Keep the C++ implementation as the readable reference.

The project's research question is still a hypothesis. A literature and prior-art review is needed before calling the method new.
