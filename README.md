# Heimdall Candidate Plotter

The Heimdall Candidate Plotter is a plotting tool designed for easy and efficient offline manual inspection of Heimdall candidates.\
Firstly, it classifies given candidates from a .cand file into several categories using user defined thresholds. See the "[Categorisation settings](https://gitlab.com/houben.ljm/heim_cand_plotter/-/edit/master/README.md#categorisation-settings)" section below for explanations of the used thresholds.

Then, these categorised candidates are shown in several [overview plots](https://gitlab.com/houben.ljm/heim_cand_plotter/-/edit/master/README.md#examples), designed to quickly obtain the desired information to discriminate true single pulses from RFI. Which overview plot is shown depends on whether the path to a single .cand file is given OR to a folder containing multiple .cand files.
 1. Path to a single .cand file is given, then only one DM versus time plot is shown (**[Figure 1](https://gitlab.com/houben.ljm/heim_cand_plotter/-/edit/master/README.md#figure-1-overview-plots-of-the-candidates-in-a-single-beam)**), displaying the candidates within the given .cand file.

      *If, in this case, other .cand files are present in the directory of the given .cand file.\
      You can use the "<" and ">" keys to switch between these files.*
 2. Path to a folder containing .cand files is given, then a separate DM versus time plot is shown for each .cand file in the directory (**[Figure 2](https://gitlab.com/houben.ljm/heim_cand_plotter/-/edit/master/README.md#figure-2-overview-plots-of-the-candidates-in-the-beams-of-a-pointing)**).\
      This makes it easier to see if candidates have occurred simultaneously in several beams and are therefore likely to be RFI.
      These candidates should be classified as RFI by the coincidencer, but depending on your data this may or may not work properly.
      Hence, this feature.

      *If, in this second case, other pointings are present in the top directory, the "<" and ">" keys can be used to switch between the
      candidates of individual pointings. However, the data must be structured [in the way described below](https://gitlab.com/houben.ljm/heim_cand_plotter/-/edit/master/README.md#assumed-data-structure) for this to work properly.*

Thirdly, if the `` --interactive`` option is given, is the shown inspection plotter, a GUI to allow interaction with the plotted candidates.
For instance, you can click on any candidate to get more detailed information about that specific candidate.
If the ``--waterfall`` option is given AND the corresponding Filterbanks are present in the pointing directory, one can also immediately inspect the RAW data of a specific selected candidate. Three options are therefore given in the plotter:
 * **"Check DM"**       - useful to check if the given DM by Heimdall is actually the value that returns a signal with the highest SNR.
 * **"Waterfall"**      - shows the dedispersed dynamic spectrum of a candidate.
 * **"Waterfall DM=0"** - shows the dispersed dynamic spectrum of a candidate (the raw data).

The plots are then saved as pngs in heim_cand_plotter/saved_plots.

These functionalities should help a user to faster determine the "realness" of candidates produced by Heimdall.
Hopefully you will find this plotter useful for your work. Happy singe pulse searching!

## Examples
#### Figure 1: Overview plots of the candidates in a single beam
![Heimdall candidate overview of a single beam.](/images/Single_beam_plot.png "Figure 1: Candidates overview of a single beam.")\
If a .cand file is given for ``cand_path`` together with the option ``--waterfall``, this is the output figure of ``show_candidates.py``.\
It is composed out of 3 sub-plots:
1.  **Time-DM plot**; shows the categorised candidates at their respective time and DM. The size of the data points shows the SNR at which the candidate was observed. Note: to not clutter this plot, is the size of these points truncated to ``SNR_THR`` set with ``--snr_thr``.
2.  **SNR-DM plot**; shows the categorised candidates at the respective true SNR and DM. SNRs > 100 are shown at an SNR of 100 to allow for better discrimination of candidates with lower SNR values. The data point sizes do not convey information.

    ***Note 1**: the colours in sub-plot 1 and 2 indicate the width of the boxcar with which a candidate was found. Which colour corresponds to which width can be seen in the colour-bar at the right of sub-plot 2.*

    ***Note 2**: the interpretation of the marker shapes is explained in the "[Categorisation settings](https://gitlab.com/houben.ljm/heim_cand_plotter/-/edit/master/README.md#categorisation-settings)" section below.*
3.  **Historgram plot**; shows how many candidates have a DM within the range of a specific DM bin. Therefore, is the full DM range linearly divided into ``NBINS`` bins, a number which can be set with the option ``--nbins``. Since the output .cand file of the ``coincidencer`` is shown here, a histogram of each beam is visible. Otherwise this sub-plot will only show the histogram plot of the candidates in the loaded beam.

Also shown are three clickable buttons (only visible due to the set option ``--waterfall`` and explained for Figure 2), a selection region that shows more detailed information about a selected candidate and which .cand file was loaded.

#### Figure 2: Overview plots of the candidates in the beams of a pointing
![Heimdall candidate overview of multiple beams.](/images/Multi_beam_plots.png "Figure 2: Candidates overview of multiple beams.")\
If for ``cand_path`` the path to a folder is given that contains the .cand files for each beam of a pointing, this is the output figure of ``show_candidates.py``.
It is composed out of multiple sub-plots similar as those in Figure 1.
1.  **Time-DM plots**; if present, the output .cand file of the ``coincidencer`` is shown in the first Time-DM sub-plot at the top of this figure. All subsequent Time-DM plots correspond to the candidates from a specific beam, indicated by the number in front of the sub-plot.
2.  **SNR-DM plot**; only visible when the ``coincidencer`` output .cand file is present. The plot has been rotated 90 degrees compared to the sub-plot in Figure 1 to align the DM values with the Time-DM plot in front of it. This eases the association of candidates in both sub-plots.
3.  **Histogram plots**; for each beam a separate histogram plot is made to show how many candidates fall within a set DM bin. Useful, to identify recurring candidates at a specific DM (for instance pulsars).

The interpretation of the colours and marker shapes are the same as for Figure 1.\
The selection region here also shows more details about a selected candidate. However, in contrast to Figure 1, is now in each Time-DM plot a selection aid displayed. The selected candidate is highlighted in the coincidencer output plots and a cross-hair is shown in the others. The cross-hair helps to guide the eye to see if, at a similar time and DM as the selected candidate, another candidate was detected in other beams. Whenever this is the case, this increases the chance of the selected candidate to be coincident RFI.\
On the right hand side of the figure one can also see which pointing's candidate files were loaded.

Like in Figure 1, three buttons are visible when the option ``--waterfall`` is given:
 * **Check DM**; if clicked, opens the new figure shown in the red box. If a real pulse was selected this plot should show a cross, or butterfly shaped pattern with the centre of the shape indicating the DM at which the de-dispersed time-series has the highest SNR. Although not clearly visible in this plot, a small "cross" can be seen at a DM of about 600 pc/cm3 (in the centre of the plot).
 * **Waterfall**; opens a new figure, as shown in the green box, with dynamic spectra of the selected candidate's Filterbank data de-dispersed to its corresponding DM.
 * **Waterfall DM=0**; same as "Waterfall", but now de-dispersed to a DM of 0 (blue box). I.e. not de-dispersed and showing the RAW Filterbank data "as is".

*Waterfall plots explained:*\
These figures are composed of multiple dynamic spectra, each centred around the time of the selected candidate and de-dispersed with the candidate's DM or DM = 0. If one goes down over the vertical axis, the dynamic spectra are increasingly more downsampled in frequency, where they are increasingly more downsampled in time towards the right of these figures. The downsampling factors are given to the side or below the dynamic spectra.\
On the top row are the time series shown obtained by summing all frequency channels in time for a given downsampling factor. On the right are the pulse's spectra given, obtained by summing al time samples in frequency for a specific downsampling factor. The MJD corresponding to time "0.0" in the dynamic spectra is given on the bottom of this figure.\
**Note**: the amplitude scale of the dynamic spectra is optimised to suppress all radiation that is dimmer than 2 times the RMS below the median and boost the visibility of radiation that is brighter than ``SNR_CUT`` times the RMS above the median. A value that can be set with the ``--snr_cut`` option.

## Getting Started
### Prerequisites
This Heimdall candidate plotter is a python 3 code and uses the modules ``matplotlib``, ``mpl_toolkits``, ``numpy`` and ``scipy``, which might not be installed with your standard python3 build.\
It further depends on the adapted ``pysigproc`` module, originally written by [Devansh Agarwal](https://github.com/devanshkv/pysigproc), to be able to import Filterbank data. This and all other dependencies are however present in this repository.

### Installing
Install python 3 with the above modules and go find single pulses!\
Though, first ensure that your data is in the correct (assumed) data structure!!!

### Assumed data structure:
In the included test folder one can see an example of the assumed data structure.\
Since the code was designed for the HTRU-North survey, data is structured based on a pointing's name and MJD at which it was observed. Though, the naming is not that important as long as your Filterbank files reside in a sub-sub-directory from the base directory like:\
```<base-dir>/<pointing name>/<MJD>/*_<beam num>_<something>.fil```

The name of your Filterbank files should have the form of: ```*_<beam num>_<something>.fil```\
In order for the code to determine which Filterbank files belong to which beam.

In the given data folder, one should have created one .cand file, per beam, that contains all the candidates produced by Heimdall for that beam. Move these combined .cand files to a separate directory.\
**Note**: it is *important* that this directory containing the combined candidate files per beam has a name that ensures it to be sorted on the bottom of a list of all the folders in your pointing directory.

Lastly, if the coincidencer is used, the coincidenced candidate file should have a name like: ```<date>-<time>_coinced_all.cand```

To make this all clear, hereby an example with a dataset of a pointing observed using three beams which resides in the base directory ```~/Test-data```.\
Per Filterbank file (i.e. per beam) a folder was created to store the short candidate files produced by Heimdall. These short .cand files where combined
in a single ``*_<beam_num>_all.cand`` file and moved to the folder ``X-all_cand``. This requered folder has been given a name such that in an alphabetical
ordering of the directories in ``~/Test-data/test_pointing/58946`` it will show up last.\
Now, one can either get the output from [Figure 1](https://gitlab.com/houben.ljm/heim_cand_plotter/-/edit/master/README.md#figure-1-overview-plots-of-the-candidates-in-a-single-beam) by setting ``cand_path`` to an individual .cand file, for instance ``~/Test-data/test_pointing/58946/X-all_cand/2015-09-21-14_45_48_02_all.cand``.
Or, get the output from [Figure 2](https://gitlab.com/houben.ljm/heim_cand_plotter/-/edit/master/README.md#figure-2-overview-plots-of-the-candidates-in-the-beams-of-a-pointing) by setting ``cand_path`` to ``~/Test-data/test_pointing/58946/X-all_cand``.
```
~/Test-data/test_pointing/58946/B0329_01_8bit.fil
~/Test-data/test_pointing/58946/B0329_01_8bit_cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_45_48_01.cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_46_06_01.cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_46_34_01.cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_47_03_01.cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_47_32_01.cand

~/Test-data/test_pointing/58946/B0329_02_8bit.fil
~/Test-data/test_pointing/58946/B0329_02_8bit_cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_45_48_02.cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_46_06_02.cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_46_34_02.cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_47_03_02.cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_47_32_02.cand

~/Test-data/test_pointing/58946/B0329_03_8bit.fil
~/Test-data/test_pointing/58946/B0329_03_8bit_cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_45_48_03.cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_46_06_03.cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_46_34_03.cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_47_03_03.cand
~/Test-data/test_pointing/58946/B0329_01_8bit_cand/2015-09-21-14_47_32_03.cand

~/Test-data/test_pointing/58946/X-all_cand
~/Test-data/test_pointing/58946/X-all_cand/2015-09-21-14_45_48_01_all.cand
~/Test-data/test_pointing/58946/X-all_cand/2015-09-21-14_45_48_02_all.cand
~/Test-data/test_pointing/58946/X-all_cand/2015-09-21-14_45_48_03_all.cand
~/Test-data/test_pointing/58946/X-all_cand/2015-09-21-14_45_48_coinced_all.cand
```

## Usage:
```
show_candidate.py [-h] [-f FIL_FILE] [-w] [-v] [--snr_cut SNR_CUT]
                  [--mem_cut MEM_CUT] [--dm_cut DM_CUT]
                  [--nbeam_cut NBEAM_CUT] [--filter_cut FILTER_CUT]
                  [--beam_mask BEAM_MASK] [--chan1 CHAN1]
                  [--chanbw CHANBW] [--tsamp TSAMP] [--nbins NBINS]
                  [--snr_thr SNR_THR] [--mps MPS] [--duration DURATION]
                  cand_file
```

### Positional arguments:
```
cand_path               Path to candidate file or folder.
```

### Optional arguments:
```
  -h, --help            show the help message
  -f FIL_PATH, --fil_file FIL_PATH
                        Path to Filterbank file or folder. (default: None)
  -w, --waterfall       Enable Filterbank plot tools. If True & 'fil_path' =
                        None the later defaults to '../'. (default: False)
  -v, --verbose         Print more operation details (default: False)
  -i, --interactive     Launches GUI (default: False)
```

### Categorisation settings:
```
  --snr_cut SNR_CUT     SNR cut below which cands are classified as hidden.
                        (default: 6.5)
  --mem_cut MEM_CUT     Members cut below which cands are classified as noise.
                        (default: 3)
  --dm_cut DM_CUT       DM cut below which cands are classified as noise.
                        (default: 1.5)
  --nbeam_cut NBEAM_CUT
                        Number of beams above which cands are classified as
                        coincident noise. (default: 3)
  --filter_cut FILTER_CUT
                        Boxcar width in units of '2**filter' above which cands
                        are classified as too fat. (default: 10)
  --beam_mask BEAM_MASK
                        Beam mask from coincidencer. (default: 32767)
  --chan1 CHAN1         Highest freq. (in MHz) of data used to classify cands
                        as too narrow. (default: 1510)
  --chanbw CHANBW       Channel bandwidth (in MHz) of data used to classify
                        cands as too narrow. (default: 0.5859375)
  --tsamp TSAMP         Sampling time (in sec) of data used to classify cands
                        as too narrow. (default: 5.4613e-05)
```

From the above categorisation setting instructions it follows that `show_candidates.py` groups Heimdall candidates into 5 categories:
1.  **hidden**; candidates with a lower SNR than `SNR_CUT` are classified as "hidden" and not shown in the Time-DM and SNR-DM plots.
2.  !["Noise" marker](/images/noise.png "") = **noise**; all candidates classified as "noise", by either having a too low DM or consist out of too few members, are shown with a grey cross.
3.  !["Coincident" marker](/images/coinc.png "") = **coincident**; if a candidate occurs in more than `NBEAM_CUT` beams or is marked as coincident by the `coincidencer`, it is shown with a grey 'plus'.
4.  !["Wrong Width" marker](/images/wrong_width.png "") = **wrong_width**; the transparent squares indicate candidates that are either too wide or too narrow. The first is determined with the threshold `--filter_cut`, the later is defined as pulses that are narrower than the dispersion smear in the highest frequency channel.
5.  !["Valid" marker](/images/valid.png "") = **valid**; all candidate not classified as "hidden", "noise", "coincident" or having a "wrong width" are classified as possible valid candidates and indicated with solid circles.

***Note**: the colours of the "wrong width" and "valid" candidates show the widths' of these candidates.*

### Histogram plot settings:
```
  --nbins NBINS         Number of bins used to make the histogram plots. (default: 30)
```

### SNR-DM plot settings:
```
  --tsamp taken from 'Categorisation settings'
```

### Time-DM plot settings:
```
  --snr_cut taken from 'Categorisation settings'.

  --mps MPS             Max Point Size of markers. (default: 30.0)
  --snr_thr SNR_THR     SNR threshold above which marker sizes are reduced to MPS.
                        (default: 25.0)
  --duration DURATION   Extend Time-DM plot from 0. to '--duration' seconds.
                        (default: None)
```

## Run a test:
First, ensure that python 3 is executed by default.\
Second, cd into the git repo which contains this README.\
Then run either one of the following two commands to launch the plotter on the provided test candidate files.\
Note, that their corresponding Filterbank files are not provided, so one can not inspect the raw data, or in other words, use the ```--waterfall``` option.

To plot a single candidate file:
```
python ./show_candidates.py -v -i test_cands/Z023.5-01.7/57286/all_cand/2015-09-21-14_45_48_01_all.cand
```

To download the png for a single candidate file without launching the interactive plotter:
```
python ./show_candidates.py -v test_cands/Z023.5-01.7/57286/all_cand/2015-09-21-14_45_48_01_all.cand
```

To plot all candidates from a pointing:
```
python ./show_candidates.py -v -i test_cands/Z023.5-01.7/57286/all_cand
```


## Notes:
 * In the event you encounter problems with the plotter, feel free to contact the author, so we can work them out together.
 * Classification and plotting based on:
   https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Scripts/trans_gen_overview.py

 * When scrolling through pointings using the "<" and ">" keys, one might get the following error when closing the plotter:
   ```
   Exception ignored in: <function WeakMethod.__new__.<locals>._cb at 0x7fd88188a050>
   Traceback (most recent call last):
    File "<...>/lib/python3.7/weakref.py", line 55, in _cb
   AttributeError: 'NoneType' object has no attribute '_alive'
   ```

   This is an unresolved bug in matplotlib as of 07.04.2020. To get rid of this message, a workaround has been suggested [here](https://github.com/matplotlib/matplotlib/issues/13033). Though, this error should not hinder the proper working of the plotter.

## Finally
### Author
*  [**Leon Houben**](https://www.linkedin.com/in/houbenljm) - Radboud University / Max Planck for Radio Astronomy

### Acknowledgements
Big thanks to Luc Builtjes who showed that Heimdall's coincidencer does not work for non-synchronised beams, which let to the development of
the multi-beam view of this plotting tool.
