# Heimdall Candidate Plotter

The Heimdall Candidate Plotter is a plotting tool designed for easy and efficient offline manual inspection of Heimdall candidates.\
First, it classifies given candidates from a .cand file into several categories using user defined thresholds. See the "[Categorisation settings](https://gitlab.com/houben.ljm/heim_cand_plotter/-/edit/master/README.md#categorisation-settings)" section below for explanations of the used thresholds.

Then, these categorised candidates are shown in several [overview plots](https://gitlab.com/houben.ljm/heim_cand_plotter/-/edit/master/README.md#examples), designed to quickly obtain the desired information to discriminate true single pulses from RFI. The path to a single .cand file is given and the candidates within the file are displayed (**[Figure 1](https://gitlab.com/houben.ljm/heim_cand_plotter/-/edit/master/README.md#figure-1-overview-plots-of-the-candidates-in-a-single-beam)**).

      *If, in this case, other .cand files are present in the directory of the given .cand file.\
      You can use the "<" and ">" keys to switch between these files.*


Third, if the `` --interactive`` option is given, is the shown inspection plotter, a GUI to allow interaction with the plotted candidates.
For instance, you can click on any candidate to get more detailed information about that specific candidate.

If the `` --interactive`` option is not given, the plots are saved as pngs in heim_cand_plotter/saved_plots.

These functionalities should help a user to faster determine the "realness" of candidates produced by Heimdall.
Hopefully you will find this plotter useful for your work. Happy singe pulse searching!

## Examples
#### Figure 1: Overview plots of the candidates in a single beam
![Heimdall candidate overview of a single beam.](/images/Single_beam_plot.png "Figure 1: Candidates overview of a single beam.")\
If a .cand file is given for ``cand_path``, this is the output figure of ``show_candidates.py``.\
It is composed out of 3 sub-plots:
1.  **Time-DM plot**; shows the categorised candidates at their respective time and DM. The size of the data points shows the SNR at which the candidate was observed. Note: to not clutter this plot, is the size of these points truncated to ``SNR_THR`` set with ``--snr_thr``.
    **Time-Beam plot**; shows the categorised candidates at the respective time and beam. The size of the data points shows the SNR at which the candidate was observed. (This plot is not shown in the figure.)
2.  **SNR-DM plot**; shows the categorised candidates at the respective true SNR and DM. SNRs > 100 are shown at an SNR of 100 to allow for better discrimination of candidates with lower SNR values. The data point sizes do not convey information.

    ***Note 1**: the colours in sub-plot 1 and 2 indicate the width of the boxcar with which a candidate was found. Which colour corresponds to which width can be seen in the colour-bar at the right of sub-plot 2.*

    ***Note 2**: the interpretation of the marker shapes is explained in the "[Categorisation settings](https://gitlab.com/houben.ljm/heim_cand_plotter/-/edit/master/README.md#categorisation-settings)" section below.*

3.  **Histogram plot**; shows how many candidates have a DM within the range of a specific DM bin. Therefore, is the full DM range linearly divided into ``NBINS`` bins, a number which can be set with the option ``--nbins``. Since the output .cand file of the ``coincidencer`` is shown here, a histogram of each beam is visible. Otherwise this sub-plot will only show the histogram plot of the candidates in the loaded beam.

The three buttons visible in this figure and Figure 1 are not shown in DSA graphs because they use a different style of cand file.

## Getting Started
### Prerequisites
This Heimdall candidate plotter is a python 3 code and uses the modules ``matplotlib``, ``mpl_toolkits``, ``numpy`` and ``scipy``, which might not be installed with your standard python3 build.\
It further depends on the adapted ``pysigproc`` module, originally written by [Devansh Agarwal](https://github.com/devanshkv/pysigproc), to be able to import Filterbank data. This and all other dependencies are however present in this repository.

### Installing
Install python 3 with the above modules and go find single pulses!\

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

To plot a single candidate file:
```
python ./show_candidates.py -v -i test_cands/heimdall.cand
```

To download the png for a single candidate file without launching the interactive plotter:
```
python ./show_candidates.py -v test_cands/heimdall.cand
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
