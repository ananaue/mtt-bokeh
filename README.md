# MTT-Bokeh

Generate graphs and IC50s from MTT raw .xlsx data in bulk.

- Supports only the "pchip" (Piecewise Cubic Hermite Interpolating Polynomial) regression currently. Functionality for "akima", "4PL", "5PL", "b-spline", "catmull_rom", and/or straight lines will be added back in the future.
- Written without the use of AI tools.
- A small set of example data is provided in `./rawdata`. It generates the following example output:

<p align="center">
  <img src="output/2026-02-08 19,27,52.019725/image.png" width="700" alt="Example output of MTT-Bokeh, showing a grid of 6 4-drug cytoxicity curves with IC50 values in their legends.">
</p>

## 1. Requirements
- Install selenium and either gecko or chromedriver, as discussed under ["install with pip" (unless using conda)](https://docs.bokeh.org/en/latest/docs/user_guide/output/export.html). If using windows, Add the folder containing the driver in Path under User Environment Variables (different process for mac/linux)
  - This is needed for PNG and SVG export.
- Install packages listed in requirements.txt (run `pip install -r requirements.txt`)

Tested using Python 3.13.5, with the latest versions of requirements.txt packages as of Feb 2, 2026.

## 2. Create template file(s).
Named either `template.xlsx` or `template#.xlsx`. For each well, a number of `[ID].pos` absorbance values are added, with (optional) `[ID].neg` values subtracted. In the case of MTT, `A630nm` is subtracted from `A562nm`. 
- Importantly, the well absorbance data must begin on row 38, column C, with respective concentration values to the left, in column B. 
  - Autodetection of where the plate data is in the template file should be added in the future, but for now one can change these properties in the definitions of `tp` and `md` in `main.py`.

- `TF` is the tag for changing template file. 
  - Ex., `+TF=1` will switch to `template1.xlsx`
  - If not added, `template.xlsx` is used by default.
- The template may not need to be updated upon missing sample(s), if they are the last on plate. 
- I have not had the chance to see if this system can be adapted for plate sizes other than 8x12.

## 3. Put properly named Excel files in `./rawdata`.
Naming format example (square brackets ommited): 
```
MM.DD.YY [Cell line ID] [Sample ID]-[Sample ID]-[Sample ID]-[Sample ID] XXhrs treatment..mb+MC=H04H06+VC=H10H12+RM1=F02F02...xlsx
```
- Sample list must be deliminated by dashes (or similar non-space, non-letter/digit potentially) and directly proceed " ##hrs" and directly follow a space (" "). With the current regex, sample IDs can be 4 digits or 3 letters (assigned to the variable `drugname`. Following "..mb" (short for MTT Bokeh) and terminated with ".." are the *filetags*.

Required *filetags* are `MC` ("mortality control") and `VC` ("vitality control").
- Ex., `"..mb+MC=H04H06+VC=H07H09.."`
- `RM` tags (e.g., `RM1`, `RM2`, `RM3`, ..) will remove known outliers.
  - Set them equal to a selection of wells (`H07H09`). If only one well, repeat the well twice (`H07H07`)).
- Andvanced Renamer and similar tools are useful for adding tags.

## 4. Edit the second section of `main.py` to change various graph config values
- Currently, only pchip reggression is known to be functional. Some changes have not been made yet for the other regressions' code due to recent script changes.

## 5. Run `main.py`
- IC50s.csv output may have an issue with encoding of "±". Also, it currently does not distinguish very well between missing data vs where an IC50 could not be calculated (such as when it is not reached in the given concentration range).

Note that the analysis does not currently correct for dead cell debris background.
- Along with a media control (here, used as `MC`) and a cells-only one (`VC`), a known 100% cytotoxicity control would be ideal. Changes would need to be made to the script to support it (potentially adding a `BC` tag (background control) and changing the calculation section)