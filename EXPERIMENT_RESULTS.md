# First experiment: precision and one-ULP nudges

This is the output from the current code using 6,000 steps and a time step of 0.01 model-time units. The full row-by-row output is in [stormline.csv](stormline.csv), and the corresponding sensitivity map is [stormline.svg](stormline.svg).

## What happened

The unperturbed `double` and `float` runs both ended with zero sites above the rain threshold. Their total rain proxies were 0.6296 and 0.7240. The total absolute difference by site was 0.3433, and their final-state RMS difference was 3.6012.

Near the initial circulation value of 8, the `float` one-ULP nudge is about `9.54e-7`; the `double` nudge is about `1.78e-15`. The CSV records the exact step used in each row.

Changing one starting value by one ULP caused a site to cross the rain threshold in three of the 20 `double` runs and three of the 20 `float` runs:

| Arithmetic | Nudged site | Site whose threshold status changed |
|---|---:|---:|
| double | 1 | 5 |
| double | 2 | 3 |
| double | 13 | 8 |
| float | 0 | 8 |
| float | 6 | 9 |
| float | 17 | 3 |

In this table, “nudged site” is where the initial value changed. The final column is where the rain-threshold result changed. The two sites do not have to match because the model couples the sites around a ring.

## How to read this result

This run shows that tiny initial changes can alter threshold outcomes in this particular toy model. It does not show that real storms are triggered by computer roundoff, and it does not establish that the method is new.

The end-state differences are already large by the end of the 60-unit run. That makes the final-state metric less useful for ranking the nudged sites at this duration. A good next experiment is to record when each paired run first separates by a chosen amount, and repeat the comparison at shorter durations.

## Reproduce it

From the project folder, run:

```powershell
stormline.exe --steps 6000 --dt 0.01 --output stormline.csv
```

The output file includes the step count and time step on every row. The result can change if the model constants, compiler options, or implementation change.
