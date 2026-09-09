# NPS-analysis-v2

Modified version of the mechano-NPS analysis code published with the JoVE
protocol paper. **[README_JOVE.md](README_JOVE.md) is still the reference** for
the citation, the calling examples, and the bundled example data and results.
This file covers only what is different here.

The two main functions were renamed, so both versions can sit on the path at once:

| JoVE original | this version |
| --- | --- |
| `mNPS_procJOVE.m` | `mNPS_procData.m` |
| `mNPS_readJOVE.m` | `mNPS_readData.m` |
| — | `mNPS_showWindow.m` (new) |

**Needs:** MATLAB R2023a or later (`xregion`), Signal Processing Toolbox, Curve
Fitting Toolbox.

---

## 1. Running it

Edit the parameters at the top of [procdata.m](procdata.m) and run it:

```matlab
filepath   = '.../young_treated_trial_1.mat';  % .mat file holding a 1xn `data` vector
ch_height  = 22.3;          % channel height [µm]
De_np      = 27.65;         % calibrated effective diameter of the node-pore segments [µm]
wC         = 10.5;          % contraction pore width [µm]
thresholds = [5e-5, 2e-4];  % [bottom, top] pulse-detection thresholds [data units]
sampleRate = 50000;         % DAQ sampling frequency [Hz]
```

A save dialog opens when processing finishes, defaulting to
`<datafile>_output.xlsx` beside the input file. Cancel it and the script stops
with `output_table` still in the workspace, so nothing needs reprocessing.

### The menu

Each candidate event is drawn in figure 42, then:

```
Is anything wrong?
  ENTER  skip this pulse
  [ 1 ]  Save window 1 (highlighted)
  [ 2 ]  Adjust window
  [ 3 ]  Auto threshold
  [ 4 ]  Adjust thresholds
  [ X ]  Stop and Output
```

* **1** saves the highlighted window (its number is in the menu line).
* **2** lets you pick another: type a number to move the highlight, then you are
  back here — nothing is saved yet, so the thresholds can still be adjusted
  before **1**. ENTER on an empty selection returns here with the highlight
  where it is. Window 1 is the default, and adjusting the thresholds (**3**,
  **4**) resets to it, since new thresholds renumber the windows.
* **3** sets top = 0.85 × the auto value, bottom = 0.12 × top, floored at the
  bottom threshold you passed in.
* **4** opens a submenu (auto or manual, top or bottom). Each change is applied
  immediately — the event is re-read and replotted. ENTER returns here, same
  event.
* **X**, or anything unrecognized, confirms then stops and returns what is saved.

The old keys still work: `P . /` save, `//` choose a window, `+ = '`
auto-threshold, `T` and `Y` thresholds.

### Figure 42

Four panels on a shared x-axis (index into the filtered, downsampled signal),
top to bottom:

| panel | contents |
| --- | --- |
| `y_diff` | first difference, with both thresholds. Blue circles = segment starts, red = ends |
| `y_LP` | filtered signal and ASLS baseline. Candidate windows numbered along the top |
| `y_rect` | rectangularized pulses |
| recovery fit | `poly1` fit to the three recovery pulses of the highlighted window |

Segment colors show each segment's role **in the highlighted window**: red =
reference (sizing), yellow = contraction, blue = recovery, grey = outside it.
The green band spans the whole window.

---

## 2. What changed

### Results

Six changes alter values in `output_table`, all on events you would have been
offered either way. Re-process rather than compare old and new tables directly.
`mNPS_readData.m` section numbers refer to its `%% SECTION` headers.

---
#### `rec_rate` could come from the wrong window · §9b

`fo` is the loop variable, so after the `for k` loop it holds the *last*
window's fit, which then gets broadcast to every row. Harmless when a read
yields a single window, the usual case; wrong when it yields two or more and you
save any but the last — including window 1, the default. It also errored
outright when the loop ran zero times, since `fo` is only ever assigned inside
it.

```matlab
% ORIGINAL
% recovery rate
calculated(:,8) = fo.p1;
```
```matlab
% MODIFIED
calculated(:,8) = out(:,9);
```

`out(k,9)` was always the correct per-window slope, so in tables you have
already saved the bug shows up as `fo_p1 ~= rec_rate`:
`nnz(T.fo_p1 ~= T.rec_rate)` counts the affected rows.

---

#### The recovery fit used the wrong times · §9

The old vector had a fabricated `0`, measured two of its points from the
**start** of the contraction (`sq_k,1`), and paired each amplitude with the
**end** index (`,2`) of the preceding segment. The new one is the start index
(`,1`) of all three recovery segments measured from the **end** of the
contraction (`sq_k,2`) — the origin `T_rec` already used. Changes `fo_p1`,
`fo_p2`, `gof_rsquare` and `rec_rate`.

```matlab
% ORIGINAL
% populate vector rT with time-points of recovery pulses
rT = [ 0, ...
    (pulse_series(rec_k_start,2)   - pulse_series(sq_k,1)), ...
    (pulse_series(rec_k_start+1,2) - pulse_series(sq_k,1)) ] ...
    ./ Fs .* N; % in ms
```
```matlab
% MODIFIED
rT = ( pulse_series(rec_k_start:rec_k_end, 1) - pulse_series(sq_k,2) )' ...
    ./ Fs .* N; % in ms
```
---
#### `dT_np` was always `NaN` · §9

Hard-coded and never replaced, even though the per-segment times on the next
line were already being computed for the velocity.

```matlab
% ORIGINAL
% average node-pore transit time in reference segments [ms]
dT_np = nan;

% node-pore transit time in each reference segment [ms] (row vector)
dT_nps = ( pulse_series(ref_k_start:ref_k_end,2) - pulse_series(ref_k_start:ref_k_end,1) )' ./Fs.*N;
```
```matlab
% MODIFIED
% dT_np = nan;

dT_nps = ( pulse_series(ref_k_start:ref_k_end,2) - pulse_series(ref_k_start:ref_k_end,1) )' ./Fs.*N;
dT_np = mean(dT_nps,2);
```
---
#### `def_diameter` used an unexplained formula · §9b

The old expression carried an unexplained `0.01` and a `pi/4`. The bracket is
the De Blois-Bean inversion, returning `d_c^3` for the volume-equivalent
sphere; modelling the deformed cell as an oblate spheroid with
`V_deform = pi/6*wC*L_deform^2` gives `L_deform = sqrt(d_c^3/wC)` — no cube
root, no `0.01`, no `pi/4`.

```matlab
% ORIGINAL
% deformed diameter (um)
calculated(:,5) = 0.01*(pi/4*wC)*(((out(:,4)./out(:,2)*De_c^2*L)./ ...
    (1+0.8*L/De_c*out(:,4)./out(:,2))).^(1/3)).^2;
```
```matlab
% MODIFIED
calculated(:,5) = sqrt( ((out(:,4)./out(:,2)*De_c^2*L)./ ...
    (1+0.8*L/De_c*out(:,4)./out(:,2))) ./ wC );
```

---
#### `wCDI` is now two columns · §9b and `mNPS_procData.m`

The arithmetic in `mNPS_readData` is **unchanged** — `wCDI_original` is the
JoVE quantity under a new name, each cell divided by its own `V_np`. What is
new is a second column dividing every cell by one flow velocity, which can only
be formed once the whole run has been picked over. See §3 for how they differ in
use.

```matlab
% ORIGINAL — one column
% wCDI (dimensionless)
calculated(:,6) = calculated(:,4)./calculated(:,3) .* calculated(:,1)./ch_height;

calc_cols = {..., 'wCDI', ...};
```
```matlab
% MODIFIED — identical formula, renamed
calculated(:,6) = calculated(:,4)./calculated(:,3) .* calculated(:,1)./ch_height;

calc_cols = {..., 'wCDI_original', ...};
```
```matlab
% ADDED — mNPS_procData.m, once output_table is built
U_flow = mean(output_table.V_np, 'omitnan'); % [mm/s]
output_table.wCDI_mean_Uflow = ...
    output_table.V_c ./ U_flow .* output_table.diameter ./ ch_height;
```
---
### Runs that used to be lost

* **The last candidate was never processed** — `while (i < size(uni_win,1))`
  dropped `uni_win(end)`, and a recording with exactly one candidate produced
  nothing at all. Now `<=`, likewise in `mNPS_cleanKim.m`.
* **A sparse recording crashed misleadingly.** `out_cols`, `out_units` and
  `rec_cat_description` were assigned inside the per-window loop, so fewer than
  6 detected segments left them undefined and MATLAB threw *"Output argument
  'rec_cat_description' is not assigned"* — from the first call, which sits
  outside the `try`/`catch`. They are constants and now set once, above the loop.
* **The contraction search could index past the end of the array**, and the
  opposite case (`find` returning `[]`) silently skipped the zeroing. Both are
  handled explicitly.
* **A typo at the main prompt ended the run**, dumping whatever was saved. It now
  asks first.
* **Threshold entry could end the run or disable thresholding.** ENTER gave `[]`,
  read as "unrecognized"; junk gave `NaN`, and `abs(x) < NaN` is never true, so
  every sample read as a pulse. Entry is now read as a string and validated.
* **Menu input is trimmed**, and ENTER matches `case ''` — the original tested
  `case []` against a string, so a trailing space fell through to stop-and-dump.
* **The completion beep could kill the run** — `sound()` ran before the table was
  returned. Now wrapped in `try`/`catch`.
* **Window selection validates its input**, and `case '1'` no longer renumbers
  its table (its "Window 1" could name a different window from the one
  highlighted and saved).

### Interface

* Numbered menus, plus a threshold submenu that applies each change immediately
  and leaves the last good result on screen if a threshold finds nothing.
* Window selection is select → highlight → confirm, staying on the same event so
  windows can be compared without re-reading.
* **Figure 43 is gone** — the recovery fit is the fourth panel of figure 42, so
  it belongs to the window on screen rather than to whatever was processed last.
  Panels are positioned directly instead of through `subplot`;
  `figwin_tighten.m` is no longer used.
* Candidate windows are numbered in the `y_LP` panel, matching the selection
  table. Segment colors show role rather than alternating by parity.
* `mNPS_readData` returns `win_range`, the `[start, end]` index range of each
  candidate window.
* **The rectangularize pass runs once, not 50 times** — it reads `y_detrend` and
  `nz_mat` and modifies neither, so 49 repeats were no-ops.
* [procdata.m](procdata.m) saves through a `uiputfile` dialog to `.xlsx`, instead
  of `save('output_table_001.mat', ...)`.

---

## 3. Reference

Inputs are not documented in README_JOVE, so they are here. Everything from
`thresholds` on is optional.

```matlab
[output_table, filtered_data, y_baseline, t_filtered, fs_filtered] = ...
    mNPS_procData(filepath, ch_height, De_np, wC, thresholds, sampleRate, ASLS_param, eventlength_filt)
```

| argument | type / units | description |
| --- | --- | --- |
| `filepath` | char | `.mat` file holding `data`, a 1×n vector of current at constant voltage [any units] |
| `ch_height` | double [µm] | channel height, from the SU-8 wafer |
| `De_np` | double [µm] | effective diameter of the node-pore segments, from calibration particles |
| `wC` | double [µm] | contraction channel width |
| `thresholds` | 1×2 double [data units] | `[bottom, top]`. Default `[1e-4, 1e-3]`; `[]` for the default |
| `sampleRate` | int [Hz] | default `50e3` |
| `ASLS_param` | struct | for `ASLS.m`. Default `lambda=1e9, p=3e-3, noise_margin=1e-4, max_iter=20` |
| `eventlength_filt` | int | filtered samples expected for one transit. Default `2000` |

Outputs are as described in README_JOVE, with the table now n×20.

```matlab
[OUT_array, empty, auto_thresh_value, column_names, column_units, rec_cat_description, ...
 filtered_data, y_baseline, t_filtered, fs_filtered, win_range] = ...
    mNPS_readData(data_vector, sampleRate, ch_height, De_np, wC, thresholds, plotflag, fitflag, ASLS_param)
```

Same geometry arguments, plus `plotflag` (draw figure 42) and `fitflag` (run the
recovery fit). `mNPS_procData` calls it three times per event: once over the
whole recording to find candidates, once on a 1-second slice with both flags off
to test for a pulse, and once with both on to display it. `auto_thresh_value` is
`[]` unless `plotflag` is true; `win_range` rows align 1:1 with `OUT_array`.

`mNPS_showWindow(window_index)` highlights one candidate window in figure 42.
The handles and fit data are stashed on the figure by `mNPS_readData`, so only
the index is needed.

### `output_table` columns

20 columns; units are also available as
`output_table.Properties.VariableUnits`.

| # | column | units | description |
| --- | --- | --- | --- |
| 1 | `start_ix` | index | start index of the event in the filtered signal |
| 2 | `I_baseline` | data units | baseline current just before the event |
| 3 | `dI_np` | data units | mean current drop in the reference (sizing) segments |
| 4 | `dI_c` | data units | current drop in the contraction segment |
| 5 | `dI_c_std` | data units | standard deviation of the contraction current drop |
| 6 | `dT_np` | ms | mean transit time across the reference segments |
| 7 | `dT_c` | ms | transit time through the contraction segment |
| 8 | `fo_p1` | — | slope of the `poly1` recovery fit |
| 9 | `fo_p2` | — | intercept of the recovery fit |
| 10 | `gof_rsquare` | — | R² of the recovery fit |
| 11 | `diameter` | µm | cell diameter |
| 12 | `strain` | — | `(diameter − wC)/diameter` |
| 13 | `V_np` | mm/s | velocity in the node-pore segments |
| 14 | `V_c` | mm/s | velocity in the contraction segment |
| 15 | `def_diameter` | µm | elongation length of the deformed cell in the contraction segment; transverse deformation is `def_diameter/diameter` |
| 16 | `wCDI_original` | — | whole-cell deformability index, **normalized by this cell's own `V_np`**; inversely related to Young's modulus |
| 17 | `rec_time` | ms | time to recover after deformation, from the cell leaving the contraction segment to entering the first segment where it had recovered. `Inf` if it never recovered in the channel |
| 18 | `rec_rate` | — | per-event recovery-fit slope (same value as `fo_p1`) |
| 19 | `rec_cat` | categorical | first segment where the cell had fully recovered: `0` instant, `1`-`2` transient, `3` never recovered in the channel |
| 20 | `wCDI_mean_Uflow` | — | the same index, **normalized by one `U_flow` = mean `V_np` over every accepted cell**. Appended by `mNPS_procData`, so it is absent from `OUT_array` |

Columns 16 and 20 are the same quantity, `V_c / U * diameter / ch_height`,
differing only in `U`. `wCDI_original` (the JoVE convention, renamed from
`wCDI`) divides each cell by its own `V_np`, so it is reproducible cell by cell.
`wCDI_mean_Uflow` divides the whole run by a single flow velocity, so its values
move if you accept a different set of events; each run prints
`U_flow = ... mm/s over N cell(s)`.

`T_rec` is in `OUT_array` but dropped from `output_table` — `rec_time` is the
same thing. "Recovered" means the pulse is back within 8% of its
pre-deformation amplitude.
