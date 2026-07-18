====================================================================
GROUND HEAT MODEL - README
====================================================================

This package contains a 1D vertical heat-transfer model for paved/urban
surfaces, the environmental forcing it runs on, and two scripts for
comparing the runs and analyzing the results.

I have set things up so that everything works from one folder with very
little to edit. Below is what each file does and how to use it.


--------------------------------------------------------------------
FILES
--------------------------------------------------------------------

environmental_forcing.mat   Atmospheric forcing the model runs on (provided).
ground_heat_model.m         Main model. Runs scenarios, writes one file per run.
compare_and_save_errors.m   Compares averaged individuals vs equivalent model,
                            makes the figures, stores error metrics.
analyze_error_metrics.m     Reads stored metrics, plots error trends.

Put all files in the same folder. The scripts find each other and the
forcing file on their own, so you do not need to set absolute paths.

You need MATLAB R2021a or newer. No extra toolboxes are required.


--------------------------------------------------------------------
RUNNING THE MODEL
--------------------------------------------------------------------

Open ground_heat_model.m. Everything you normally change is near the top:

  model_dt_seconds   the model time step, in seconds. This is the only
                     timing knob. The native step of the forcing is read
                     from the file and the inputs are resampled for you.
                     The total duration also comes from the forcing.

  dz_idx             which vertical resolution to use (an index into the
                     matrix_dx list of grid spacings, in metres).

  depth_list         the domain depth(s), in metres. One value or a range.

  A_rows             which row(s) of the areal-ratio matrix A to use.
                     Each row is one mix of the six base materials.

  mat_idx            which case to run. Indices 1-6 are the individual
                     materials; indices 7-16 are the combination models.

These are nested loops, so the model runs every combination you ask for.
Press Run. Spin-up progress prints to the console.

Note: the solver is explicit, so if you make the grid finer you usually
need a smaller time step. A warning prints if the step is too large.


--------------------------------------------------------------------
OUTPUT
--------------------------------------------------------------------

Results are written to:

  Results/<today's date>/<scenario>.mat

The filename encodes the run, for example:

  m03_d0.50_dz0.010.mat                     (single material 3)
  m12_a10-20-30-10-20-10_d0.50_dz0.010.mat  (combination 12, that ratio)

Each file holds the temperature profile, the ground heat flux, the depth,
the time vector (seconds), and a meta struct with all the run settings.


--------------------------------------------------------------------
MAKING THE FIGURES
--------------------------------------------------------------------

1) compare_and_save_errors.m
   Set results_dir to your dated results folder, and set the ratio row,
   resolution index, depth and the equivalent-model index to match runs
   you already produced. It needs the six individual runs plus one
   combination run for that scenario. Running it makes three figures
   (surface temperature and flux, the depth-time profiles, and the
   individual material profiles) and adds an entry to
   Results/error_metrics.mat. Run it once per scenario you want analyzed.

2) analyze_error_metrics.m
   Reads error_metrics.mat after you have built up several entries.
   It plots error against resolution (pick one areal ratio or average
   over all) and error against material-property variance (pick one
   resolution). The trends are only meaningful once you have runs across
   several resolutions and/or several mixes.


--------------------------------------------------------------------
REGENERATING THE FORCING
--------------------------------------------------------------------

You only need this if the input data changes. Edit the paths and column
map at the top of save_environmental_forcing.m and run it once. It writes
a fresh environmental_forcing.mat in the same folder, and the rest of the
pipeline adapts to the new duration and time step automatically.


====================================================================
LICENSE AND USE REQUIREMENTS
====================================================================

This package contains both model/software files and data/materials.
Different licenses apply to these components.

MODEL AND SOFTWARE
All model code, MATLAB scripts and functions are licensed under the
PolyForm Noncommercial License 1.0.0. They may be used, copied, modified
and redistributed for noncommercial purposes only (scientific research,
education, personal study, testing and related noncommercial academic
use). This restriction also applies to any modified, extended, translated
or derivative version. No commercial use is permitted without prior
written permission from the copyright holders.

DATA AND MATERIALS
All data files, model outputs, input files, tables and figures are
licensed under the Creative Commons Attribution-NonCommercial-ShareAlike
4.0 International License (CC BY-NC-SA 4.0). They may be shared and adapted
for noncommercial purposes with appropriate attribution, and any adapted
material must be shared under the same terms. This also applies to any
modified or derivative version. No commercial use is permitted without
prior written permission from the copyright holders.

COMMERCIAL USE
Commercial use is not permitted for any part of this package unless prior
written permission has been granted by the copyright holders. Permission
for noncommercial use does not imply permission for commercial use.

THIRD-PARTY MATERIALS
Any third-party data or software included or referenced here remains
subject to its own original license and terms of use.

CONTACT
License inquiries: 
Elie Bou-Zeid
ebouzeid@princeton.edu

Technical support:
Erfan Hosseini
eh9561@princeton.edu

If there is any inconsistency between this summary and the full license
texts, the full license texts govern.