#!/usr/bin/env python

###############################################################################
#
# This code is intended to be used to quickly inspect heimdall candidates, by
# not only showing the candidates in a DM verus time plot but also dynamic
# spectra of an individually selected candidate.
#
# NOTE: for the code to properly work a fixed folder structure is assumed!!!
#       "<path to data>/<pointing-name>/<MJD>/<*_cand folders>"
#
###############################################################################

import time
import glob
import os,sys
import argparse
import subprocess
import numpy as np
import fil_cand_cut as fcc
import plotfunctions as pf
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button
from Cand_Classifier import load_candidates as LC
from Cand_Classifier import Cand_Classifier as CC
from mpl_toolkits.axes_grid1 import make_axes_locatable


class CandBrowser(object):
    """
    Functionality taken fom:
    https://matplotlib.org/examples/event_handling/data_browser.html
    """

    def __init__(self,Pcand,Pfil=None,waterfall=False,verbose=False,interactive=False):
        """
        Constructor to build plot and take care of
        dynamic repopulation of plots.
        """
        # Static settings:
        self.nbins    = 30
        self.tsamp    = 54.613e-6
        self.snr_cut  = 6.5
        self.snr_thr  = 25.
        self.mps      = 30.
        self.duration = None
        self.classify_args = []

        # Dynamic settings:
        self.Pcand = Pcand
        self.Pfil  = Pfil
        self.waterfall = waterfall
        self.verbose   = verbose
        self.interactive = interactive

    def start_browser(self):
        self.load_data()
        self.get_cands_shift()
        self.setup_figure()
        self.populate_axis()
        self.activate_picking()

        plt.savefig(os.path.join("saved_plots", os.path.basename(self.Fcands[0]) + ".png"))
        if(self.interactive):
            plt.show()

    def activate_picking(self):
        # Activate interactive picking and buttonsm
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.canvas.mpl_connect('key_press_event', self.onpress)
        if self.waterfall:
            self.dm_button.on_clicked(self.onclick_check_dm)
            self.wf_button.on_clicked(self.onclick_check_downsamp_plots)
            self.wf0_button.on_clicked(self.onclick_check_downsamp_dm0_plots)

    def load_data(self):
        ### Get data path
        self.Pcand = os.path.abspath(self.Pcand)

        # If a directory path is given enable pointing scrolling
        if os.path.isdir(self.Pcand):
            print("Here:",self.Pcand)
            #self.Pcand = 'test_cands/Z023.5-01.7/57286/all_cand'
            #self.Pcand = 'C:\\Users\\Mirabai\ Smoot\\Desktop\\DSA\\heim_cand_plotter\\test_cands\\Z023.5-01.7\\57286\\all_cand'
            #NOTE: this assumes a dir structure of <pointing name>/<MJD>/<cand dirs>
            # paths = sorted(subprocess.check_output(['find', self.Pcand, \
            #                '-mindepth', '1', '-type', 'd'], encoding='ascii').split())
            #paths = os.listdir(self.Pcand)
            paths = glob.glob(self.Pcand + '\*')
            print("len(paths):",len(paths))

            if len(paths) != 0:
                self.Pcand = paths[-1]

            Ppointing  = self.Pcand.split(self.Pcand.split('/')[-1])[0]
            Ppointing  = os.path.abspath(Ppointing)
            print("Ppointing:",Ppointing)
            Ppointings = glob.glob(Ppointing+'\..\*') #I took off one \..
            # Ppointings = sorted(subprocess.check_output(['find', Ppointing+'/../..', \
            #                '-mindepth', '2', '-maxdepth', '2', '-type', 'd'], encoding='ascii').split())
            self.Ppointings = [os.path.abspath(path) for path in Ppointings]
            self.Pointscroll = True
            self.Fcandsind = self.Ppointings.index(Ppointing)

        # If a file path is given enable beam scrolling
        elif os.path.isfile(self.Pcand):
            self.Fcands = self.Pcand
            self.Pcand  = self.Pcand.split(self.Pcand.split('\\')[-1])[0]
            self.Pointscroll = False

        # Get candidate and filterbank files
        self.get_files()

    def get_files(self):
        if not self.Pointscroll: #if given single file
            Fcand = self.Fcands
            self.Fcands = [self.Fcands]
        else:
            self.Fcands = glob.glob(self.Pcand + '\..\*')

        ### Get candidate files
        # self.Fcands = sorted(subprocess.check_output(['find', self.Pcand, \
        #                            '-maxdepth', '1', '-type', 'f',\
        #                            '-name', str("*.cand")], encoding='ascii').split())
        #self.Fcands = [self.Fcands] #if single file

        if len(self.Fcands) == 0:
            raise ValueError("Could not find .cand files in %s" % self.Pcand)
        else:
            self.Fcoincedcand = False
            self.NFcands = len(self.Fcands)

            #NOTE: this is for distinguishing between coincidenced and noncoincidenced files
            #NOTE: this assumes a file name of <date>-<time>_<beam number>_<something>.cand
            if self.Pointscroll:
                self.beam_num = [x.split('_')[-2] for x in self.Fcands]
                try:
                    #NOTE: this assumes the coincidenced candidate file has a file name of <date>-<time>_coinced_<something>.cand
                    coinced_idx = self.beam_num.index('coinced')
                except ValueError:
                    print("\n    WARNING: coincidencer's .cand file could not be found!\n")
                    self.Fcoincedcand = False
                    self.NFcands = len(self.Fcands)
                else:
                    self.Fcoincedcand = True
                    Fcoinced = self.Fcands[coinced_idx]
                    del self.Fcands[coinced_idx]
                    del self.beam_num[coinced_idx]
                    self.NFcands = len(self.Fcands)
                    self.Fcands.insert(0,Fcoinced)
                    self.beam_num.insert(0,'coinced')

        if not self.Pointscroll:
            self.Fcandsind = self.Fcands.index(Fcand)

        ### Get corresponding filterbank files
        if self.waterfall == True:
            if self.Pfil is None:
                #NOTE: this assumes the filterbank files to be located in the MJD folder along with all candidate dirs
                Pfil = os.path.abspath(self.Pcand+'/../')
            elif os.path.isfile(self.Pfil):
                Pfil  = self.Pfil.split(self.Pfil.split('/')[-1])[0]
            else:
                Pfil = os.path.abspath(self.Pfil)

            self.Ffils = sorted(subprocess.check_output(['find', Pfil, \
                                            '-maxdepth', '1', '-type', 'f',\
                                            '-name', str("*.fil")], encoding='ascii').split())
            if len(self.Ffils) == 0:
                raise ValueError("Could not find .fil files in %s" % self.Pfil)

    def get_cands_shift(self):
        classifier = CC(*self.classify_args)
        self.cands     = []
        self.multibeam = []
        self.categories = []
        for Fcand in self.Fcands:
            cs,mb = LC(Fcand,self.verbose)
            self.cands.append(cs)
            self.multibeam.append(mb)
            self.categories.append(classifier.categories(cs,mb,self.verbose))

    def setup_figure(self):
        if self.Pointscroll:
            self.fig = plt.figure(figsize=(8.3,11.7)) #A4 portrait
            #https://matplotlib.org/api/_as_gen/matplotlib.pyplot.subplots_adjust.html
            self.fig.subplots_adjust(bottom=0.05, top=0.95, hspace=0.1)
            self.fig_multi_cands()
        else:
            self.fig = plt.figure(figsize=(11.7,8.3)) #A4 landscape
            self.fig_single_cand()


    def fig_single_cand(self):
        #Build cand info plots
        #gs = gridspec.GridSpec(8, 3, height_ratios=[], width_ratios=[2, 2, 1])
        gs = self.fig.add_gridspec(nrows=16, ncols=3, width_ratios=[2, 2, 1])
        self.num_dm_ax  = [self.fig.add_subplot(gs[:8,0])]
        self.snr_dm_ax  = [self.fig.add_subplot(gs[:8,1])]
        divider = make_axes_locatable(self.snr_dm_ax[-1])
        self.cbar_ax = divider.append_axes("right", size="5%", pad=0.1)
        self.dm_t_ax    = [self.fig.add_subplot(gs[8:,:2])]

        #Build text and buttons
        self.text_ax        = self.fig.add_subplot(gs[:-4,2])
        if self.waterfall:
            self.dm_box     = self.fig.add_subplot(gs[-3,2])
            self.dm_button  = Button(self.dm_box, 'Check DM')
            self.wf_box     = self.fig.add_subplot(gs[-2,2])
            self.wf_button  = Button(self.wf_box, 'Waterfall')
            self.wf0_box    = self.fig.add_subplot(gs[-1,2])
            self.wf0_button = Button(self.wf0_box, 'Waterfall DM=0')

    def fig_multi_cands(self):
        gs = self.fig.add_gridspec(ncols=9, nrows=np.clip(self.NFcands+2,9,15))

        ### Build cand info plots
        self.dm_t_ax   = []
        self.num_dm_ax = []
        self.snr_dm_ax = []

        # Build dm_t axis for coincidenced beams
        if self.Fcoincedcand:
            self.dm_t_ax.append(self.fig.add_subplot(gs[:2,:6]))
            self.snr_dm_ax.append(self.fig.add_subplot(gs[:2,6:], sharey=self.dm_t_ax[0]))
            divider = make_axes_locatable(self.snr_dm_ax[-1])
            self.cbar_ax = divider.append_axes("right", size="5%", pad=0.1)
            self.num_dm_ax.append(0)

        # Build dm_t axis for all beams
            for idx in range(2,self.NFcands+2):
                self.dm_t_ax.append(self.fig.add_subplot(gs[idx,:6], sharex=self.dm_t_ax[0], sharey=self.dm_t_ax[0]))
                self.num_dm_ax.append(self.fig.add_subplot(gs[idx,6], sharey=self.dm_t_ax[0]))

        else:
            self.dm_t_ax.append(self.fig.add_subplot(gs[2,:6]))
            self.num_dm_ax.append(self.fig.add_subplot(gs[2,6], sharey=self.dm_t_ax[0]))

            for idx in range(3,self.NFcands+2):
                self.dm_t_ax.append(self.fig.add_subplot(gs[idx,:6], sharex=self.dm_t_ax[0], sharey=self.dm_t_ax[0]))
                self.num_dm_ax.append(self.fig.add_subplot(gs[idx,6], sharey=self.dm_t_ax[0]))

        ### Build text and buttons
        self.text_ax        = self.fig.add_subplot(gs[2:-3,7:])
        if self.waterfall:
            self.dm_box     = self.fig.add_subplot(gs[-3,7:])
            self.dm_button  = Button(self.dm_box, 'Check DM')
            self.wf_box     = self.fig.add_subplot(gs[-2,7:])
            self.wf_button  = Button(self.wf_box, 'Waterfall')
            self.wf0_box    = self.fig.add_subplot(gs[-1,7:])
            self.wf0_button = Button(self.wf0_box, 'Waterfall DM=0')

    def populate_axis(self):
        #if plotting a directory instead of a specific file
        if self.Pointscroll:
            axrange = True
            for axnum,ax in enumerate(self.dm_t_ax):
                if not all(len(x)==0 for x in self.categories[axnum].values()):
                    pf.TimeDMPlot(self.dm_t_ax[axnum], self.categories[axnum], self.duration, self.snr_cut, \
                                  self.snr_thr, self.mps, self.multibeam[axnum], axrange, False)

                    if (axnum == 0) and self.Fcoincedcand:
                        pf.DMSNRPlot_multibeams(self.snr_dm_ax[axnum], self.cbar_ax, self.categories[axnum], self.tsamp)
                    else:
                        pf.DMHistPlot_multibeams(self.num_dm_ax[axnum], self.cands[axnum], self.nbins)
                else:
                    self.num_dm_ax[axnum].tick_params(axis='y', labelleft=False)
                    self.num_dm_ax[axnum].tick_params(axis='x', bottom=False, labelbottom=False)

                if (axnum == 0) and self.Fcoincedcand:
                    ax.set_ylabel('$\\rm DM\;(pc\;cm^{-3})$', size=12)
                else:
                    ax.set_ylabel(int(self.beam_num[axnum])-1, size=12,rotation=0)

                if axnum == 0:
                    axrange = False

                if ax == self.dm_t_ax[-1]:
                    ax.set_xlabel('$\\rm Time\; (sec)$', size=12)
                    ax.tick_params(axis='x', labelbottom=True)

            #https://jakevdp.github.io/PythonDataScienceHandbook/04.10-customizing-ticks.html
            self.dm_t_ax[0].yaxis.set_major_locator(plt.LogLocator())
            self.dm_t_ax[0].minorticks_on()
        else:
            if len(self.cands[self.Fcandsind]) == 0:
                self.num_dm_ax[0].cla()
                self.snr_dm_ax[0].cla()
                self.dm_t_ax[0].cla()
            else:
                pf.DMHistPlot(self.num_dm_ax[0], self.cands[self.Fcandsind], self.nbins, self.multibeam[self.Fcandsind])
                pf.DMSNRPlot(self.snr_dm_ax[0], self.cbar_ax, self.categories[self.Fcandsind],self.tsamp)
                pf.TimeDMPlot(self.dm_t_ax[0], self.categories[self.Fcandsind], self.duration, self.snr_cut, \
                             self.snr_thr, self.mps, self.multibeam[self.Fcandsind])

        # With selection highlights
        if self.Pointscroll:
            self.selected_dmt_hlines = []
            self.selected_dmt_vlines = []

            if self.Fcoincedcand:
                min_range = 1
                select_circ = True
            else:
                min_range = 0
                select_circ = False

            for axnum in range(min_range,len(self.dm_t_ax)):
                self.selected_dmt_hlines.append(self.dm_t_ax[axnum].axhline(color='r', linewidth=0.3, zorder=0, visible=False))
                self.selected_dmt_vlines.append(self.dm_t_ax[axnum].axvline(color='r', linewidth=0.3, zorder=0, visible=False))

        if not self.Pointscroll or select_circ:
            self.selected_snrdm, = self.snr_dm_ax[0].plot([0],[0], 'o', ms=12, alpha=0.4, \
                                                     color='yellow', visible=False)
            self.selected_dmt,   = self.dm_t_ax[0].plot([0],[0], 'o', ms=12, alpha=0.4, \
                                                     color='yellow', visible=False)

        # With informative text
        initial_text = 'Click on a candidate\nto get its parameters'
        self.text_ax.cla()
        self.text_ax.axis("off")
        self.text = self.text_ax.text(0.,0.,initial_text, \
                          size=10, transform=self.text_ax.transAxes)
        if self.Pointscroll:
            pointing = self.Ppointings[self.Fcandsind].split('\\')[-2]+'/'+self.Ppointings[self.Fcandsind].split('/')[-1]
            self.text_ax.text(1.1,0.55,"Loaded pointing:\n%s" \
                % pointing, fontsize=14, rotation=-90, transform=self.text_ax.transAxes)
        else:
            self.text_ax.text(1.1,0.25,"Loaded candidate file:\n%s" \
                % (self.Fcands[self.Fcandsind].split('/')[-1]), fontsize=14, \
                rotation=-90, transform=self.text_ax.transAxes)

    def gentextinfo(self):
        infostr  = 'Selected candidate:\n\n'
        infostr += 'SNR   = %.3f\n' % self.cands[self.cand_idx[0]]['snr'][self.cand_idx[1]]
        infostr += 'DM    = %.3f ($pc\,cm^{-3}$)\n' % self.cands[self.cand_idx[0]]['dm'][self.cand_idx[1]]
        infostr += 'Time  = %.3f (sec)\n' % self.cands[self.cand_idx[0]]['time'][self.cand_idx[1]]
        infostr += 'Width = %.3f (ms)\n' % (2.**self.cands[self.cand_idx[0]]['filter'][self.cand_idx[1]] * self.tsamp*1000.)
        #infostr += 'Mem  = %d\n' % self.cands[self.cand_idx[0]]['members'][self.cand_idx[1]]
        if self.multibeam[self.cand_idx[0]]:
            infostr += "Beam = %d" % self.cands[self.cand_idx[0]]['beam'][self.cand_idx[1]]
        else:
            pass
            #infostr += "Beam = %d" % (int(self.beam_num[self.cand_idx[0]]) - 1) #Correct for zero based beam numbering

        return infostr

    def remove_bandpass(self,data,indep=True):
        """Subtract the median from each channel,
           and divide by global std deviation (if indep==False), or
           divide by std deviation of each row (if indep==True).

           Input:
               data:  Numpy array of data that needs its bandpass removed.
               indep: Boolean. If True, scale each row independantly (Default: False).

           Output:
               scaled_date: The scaled data.

           Note: taken from PRESTO's spectra.py
        """
        if not indep:
            std = data.std()
        for ii in range(data.shape[0]):
            chan = data[ii,:]
            median = np.median(chan)
            if indep:
                std = np.array([chan.std()])
                std = np.where(std==0,1,std)[0]
            chan[:] = (chan-median)/std
        return data

    def get_Fil_Cut(self):
        #Manualy garbage collect to free memory
        #gc.collect()

        if self.multibeam[self.cand_idx[0]]:
            beamN = self.cands[self.cand_idx[0]]['beam'][self.cand_idx[1]]
        else:
            beamN = int(self.beam_num[self.cand_idx[0]]) - 1 #Correct for zero based beam numbering

        for fil in self.Ffils:
            #NOTE: this assumes a fixed naming schema of: *_<beam num>_<something>.fil !!!
            if int(fil.split('_')[-2]) == beamN:
                FC = fcc.Fil_Cand(fp=fil, dm=self.cands[self.cand_idx[0]]['dm'][self.cand_idx[1]], \
                                  tcand=self.cands[self.cand_idx[0]]['time'][self.cand_idx[1]],    \
                                  _filter=self.cands[self.cand_idx[0]]['filter'][self.cand_idx[1]],
                                  verbose=self.verbose)
        try:
            FC
        except NameError:
            FC = None
            if self.verbose:
                print("\n    WARNING: No fil file for beam %d was found!\n" % beamN)

        return FC

    def onclick_check_dm(self, event):
        print("Making DM-Time plot...")
        FC  = self.get_Fil_Cut()

        if FC is None:
            return

        FC.get_chunk()
        FC.fp.close()
        FC.rebin2D(1,2) # Downsample in time for speed-up
        dmt = FC.dmtime(dmsteps=101)

        #plt.ion() # Turn interactive plotting on
        fig,ax = plt.subplots()
        ax.imshow(dmt[2], aspect='auto', \
                  extent=(dmt[0].min(),dmt[0].max(),dmt[1].min(),dmt[1].max()))
        ax.set_xlabel('Time in sec from MJD$\;%0.6f$' % FC.MJDcand, size=12)
        ax.set_ylabel('$\\rm DM\;(pc\;cm^{-3})$', size=12)
        plt.title('Bow-tie plot around DM = %0.4f' % self.cands[self.cand_idx[0]]['dm'][self.cand_idx[1]])
        #plt.tight_layout()

        del FC # Free memory
        fig.show()

    def onclick_check_downsamp_plots(self, event):
        print("Making downsamp_plots with dm=%f" % self.cands[self.cand_idx[0]]['dm'][self.cand_idx[1]])
        self.make_downsamp_plots(dm=self.cands[self.cand_idx[0]]['dm'][self.cand_idx[1]])

    def onclick_check_downsamp_dm0_plots(self, event):
        print("Making downsamp_plots with dm=0")
        self.make_downsamp_plots()

    def make_downsamp_plots(self, dm=0., max_fac=128):
        FC  = self.get_Fil_Cut()

        if FC is None:
            return

        FC.get_chunk()
        FC.fp.close()
        dedisp_data = FC.dedisperse(dms=dm, cut=False)
        burst_edge  = FC.dispersion_delay()*0.4
        FC.data = None #To save memory
        #Determine downsamp limits
        Dlim = []
        pows  = 2**np.arange(15)
        ipows = np.where(pows >= max_fac)[0][0]
        for shape in dedisp_data[2].shape:
            res   = shape/pows
            ires  = np.where(res <= 64)[0][0] #Hardcoded limit of 64 data points!
            Dlim.append(np.min((ipows,ires)))

        #Build figure
        size_ratio = 2.
        hr = np.ones(Dlim[0]+2)*size_ratio
        hr[0] = 1
        wr = np.ones(Dlim[1]+2)*size_ratio
        wr[-1] = 1
        #plt.ion() # Turn interactive plotting on
        fig = plt.figure(figsize=(1.5*np.sum(wr/size_ratio),1.5*np.sum(hr/size_ratio)))
        gs = gridspec.GridSpec(Dlim[0]+2, Dlim[1]+2, height_ratios=hr, width_ratios=wr)
        for row in range(1,Dlim[0]+2):
            for column in range(Dlim[1]+1):
                if (row == 1) and (column ==0):
                    ax0 = fig.add_subplot(gs[row,column])
                    ax = ax0
                else:
                    ax = fig.add_subplot(gs[row,column], sharex=ax0, sharey=ax0)

                rebined_data = FC.rebin2D(pows[row-1],pows[column],dedisp_data)
                rebined_data = (rebined_data[0],rebined_data[1],self.remove_bandpass(rebined_data[2],False))
                vmin = np.median(rebined_data[2])-2.*np.std(rebined_data[2])
                vmax = np.median(rebined_data[2])+self.snr_cut*np.std(rebined_data[2])
                ax.imshow(rebined_data[2],aspect='auto', \
                  extent=(rebined_data[0].min(),rebined_data[0].max(), \
                          rebined_data[1].min(),rebined_data[1].max()), \
                         vmin=vmin, vmax=vmax)

                # Set downsamp factors
                if ax.is_last_row():
                    ax.set_xlabel('fac = %dx' % pows[column])
                if ax.is_first_col():
                    ax.set_ylabel('fac = %dx' % pows[row-1])

                # Build time series
                if row == 1:
                    rebined_ts   = np.sum(rebined_data[2],axis=0)
                    bin_edge     = np.where(rebined_data[0] > -burst_edge)[0][0]
                    noise_ts     = np.concatenate((rebined_ts[0:bin_edge], \
                                                   rebined_ts[-bin_edge:-1]))
                    noise_median = np.median(noise_ts)
                    noise_std    = noise_ts.std()

                    if column == 0:
                        ts0 = fig.add_subplot(gs[0,column], sharex=ax0)
                        ts = ts0
                    else:
                        ts = fig.add_subplot(gs[0,column], sharex=ax0, sharey=ts0)
                    ts.plot(rebined_data[0],(rebined_ts-noise_median)/noise_std)
                    ts.set_xlim(rebined_data[0].min(),rebined_data[0].max())
                    ts.set_ylim(-6.,60.)
                    plt.setp(ts.get_xticklabels(), visible=False)

                # Build spectra
                if column == Dlim[1]:
                    if row == 1:
                        spec1 = fig.add_subplot(gs[row,Dlim[1]+1], sharey=ax0)
                        spec = spec1
                    else:
                        spec = fig.add_subplot(gs[row,Dlim[1]+1], sharex=spec1, sharey=ax0)
                    spec.plot(np.sum(rebined_data[2],axis=1),rebined_data[1])
                    spec.set_ylim(rebined_data[1].min(),rebined_data[1].max())
                    plt.setp(spec.get_yticklabels(), visible=False)

        # Remove intra-subplot labels
        # https://github.com/matplotlib/matplotlib/issues/11036
        all_axes = fig.get_axes()
        for ax in all_axes:
            if not ax.is_first_col():
                plt.setp(ax.get_yticklabels(), visible=False)
            else:
                ax.tick_params(axis="y", labelsize=8)
            if not ax.is_last_row():
                plt.setp(ax.get_xticklabels(), visible=False)
            else:
                ax.tick_params(axis="x", labelsize=8, rotation=45)

        # Set common axes labels
        plt.subplots_adjust(top=0.95, bottom=0.15, left=0.15, right=0.95)
        fig.text(0.5, 0.03, 'Time in sec from MJD %f' % FC.MJDcand, \
                 va='center', ha='center', fontsize=14)
        fig.text(0.03, 0.5, 'Frequency in MHz', \
                 va='center', ha='center', rotation='vertical', fontsize=14)

        del FC # Free memory
        fig.show()

    def onpress(self, event):
        if event.key not in (',', '.'):
            return True
        if event.key == '.':
            inc = 1
        else:
            inc = -1

        self.Fcandsind += inc
        if self.Pointscroll:
            self.Fcandsind = np.clip(self.Fcandsind, 0, len(self.Ppointings) - 1)
            self.Pcand = sorted(subprocess.check_output(['find', self.Ppointings[self.Fcandsind], \
                                '-mindepth', '1', '-type', 'd'], encoding='ascii').split())[-1]
        else:
            self.Fcandsind = np.clip(self.Fcandsind, 0, len(self.Fcands) - 1)
        self.reload()

    def onpick(self, event):
        N = len(event.ind)
        if not N:
            return True

        # Check to which subplot event belongs
        try:
            cand_id = self.dm_t_ax.index(event.mouseevent.inaxes)
        except ValueError:
            try:
                cand_id = self.snr_dm_ax.index(event.mouseevent.inaxes)
            except:
                return True
            else:
                if self.Pointscroll:
                    xis = 'snr'
                    yis = 'dm'
                else:
                    xis = 'dm'
                    yis = 'snr'
        except:
            return True
        else:
            xis = 'time'
            yis = 'dm'

        # the click locations
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata

        if not self.Pointscroll:
            cand_id = self.Fcandsind

        # Find nearest candidate in plot
        distances = np.hypot(x - self.cands[cand_id][xis], y - self.cands[cand_id][yis])
        indmin = distances.argmin()
        self.cand_idx = (cand_id, indmin)
        self.update()

    def reload(self):
        if self.Pointscroll:
            old_fig = self.fig
            self.get_files()
            self.get_cands_shift()
            self.setup_figure()
            self.populate_axis()
            self.activate_picking()
            self.fig.show()
            plt.close(old_fig)
        else:
            self.populate_axis()
            self.fig.canvas.draw()

    def update(self):
        # Set selection highlights
        if self.Fcoincedcand or not self.Pointscroll:
            self.selected_snrdm.set_visible(True)
            self.selected_dmt.set_visible(True)
            self.selected_dmt.set_data(self.cands[self.cand_idx[0]]['time'][self.cand_idx[1]], \
                                       self.cands[self.cand_idx[0]]['dm'][self.cand_idx[1]])
            if self.Pointscroll:
                self.selected_snrdm.set_data(self.cands[self.cand_idx[0]]['snr'][self.cand_idx[1]], \
                                             self.cands[self.cand_idx[0]]['dm'][self.cand_idx[1]])
            else:
                self.selected_snrdm.set_data(self.cands[self.cand_idx[0]]['dm'][self.cand_idx[1]], \
                                             self.cands[self.cand_idx[0]]['snr'][self.cand_idx[1]])
            clipped_snr = np.clip(self.cands[self.cand_idx[0]]['snr'][self.cand_idx[1]],0.,self.snr_thr)
            norm_snr    = (clipped_snr-self.snr_cut)/(self.snr_thr-self.snr_cut)
            self.selected_snrdm.set_markersize(np.sqrt(20)*2.)
            self.selected_dmt.set_markersize(self.mps*norm_snr*1.5)

        # Draw selection lines in multibeam plots
        if self.Pointscroll:
            for ax_id in range(0,len(self.selected_dmt_hlines)):
                self.selected_dmt_hlines[ax_id].set_visible(True)
                self.selected_dmt_hlines[ax_id].set_ydata(self.cands[self.cand_idx[0]]['dm'][self.cand_idx[1]])
                self.selected_dmt_vlines[ax_id].set_visible(True)
                self.selected_dmt_vlines[ax_id].set_xdata(self.cands[self.cand_idx[0]]['time'][self.cand_idx[1]])

        # Set informative text
        self.text.set_text(self.gentextinfo())
        self.fig.canvas.draw()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='show_candidate.py', \
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('cand_path', \
            help="Path to candidate file or folder.")
    parser.add_argument('-f','--fil_path', type=str, default=None, \
            help="Path to filterbank file or folder.")
    parser.add_argument('-w','--waterfall', action='store_true', default=False, \
            help="Enable filterbank plot tools. If True & 'fil_path' = None the later defaults to '../'.")
    parser.add_argument('-v', '--verbose', action='store_true', default=False, \
            help="Print more operation details")
    parser.add_argument('-i', '--interactive', action='store_true', default=False, \
            help="Launches GUI")

    # Categorisation settings
    classify = parser.add_argument_group('Categorisation settings')
    classify.add_argument('--snr_cut', type=float, default=6.5, \
            help="SNR cut below which cands are classified as hidden.")
    classify.add_argument('--mem_cut', type=int, default=3, \
            help="Members cut below which cands are classified as noise.")
    classify.add_argument('--dm_cut', type=float, default=1.5, \
            help="DM cut below which cands are classified as noise.")
    classify.add_argument('--nbeam_cut', type=int, default=3, \
            help="Number of beams above which cands are classified as coincident noise.")
    classify.add_argument('--filter_cut', type=int, default=10, \
            help="Boxcar width in units of '2**filter' above which cands are classified as too fat.")
    classify.add_argument('--beam_mask', type=int, default=(1<<15)-1, \
            help="Beam mask from coincidencer.")
    classify.add_argument('--chan1', type=float, default=1510, \
            help="Highest freq. (in MHz) of data used to classify cands as too narrow.")
    classify.add_argument('--chanbw', type=float, default=0.5859375, \
            help="Channel bandwidth (in MHz) of data used to classify cands as too narrow.")
    classify.add_argument('--tsamp', type=float, default=54.613e-6, \
            help="Sampling time (in sec) of data used to classify cands as too narrow.")

    # Histogram settings
    hist = parser.add_argument_group('Histogram settings')
    hist.add_argument('--nbins', type=int, default=30, \
            help="Number of bins used to make histogram.")

    # DM-SNR settings
    dm_snr = parser.add_argument_group('DM-SNR plot settings', \
            "--tsamp taken from 'Categorisation settings'")

    # Time-DM settings
    time_dm = parser.add_argument_group('Time-DM plot settings', \
            "--snr_cut taken from 'Categorisation settings'.")
    time_dm.add_argument('--mps', type=float, default=30., \
            help="Max Point Size of markers.")
    time_dm.add_argument('--snr_thr', type=float, default=25., \
            help="SNR threshold above which markersizes are reduced to MPS.")
    time_dm.add_argument('--duration', default=None, \
            help="Extend Time-DM plot from 0. to '--duration' seconds.")

    args = parser.parse_args()

    CB = CandBrowser(args.cand_path,args.fil_path,args.waterfall,args.verbose,args.interactive)
    CB.nbins    = args.nbins
    CB.tsamp    = args.tsamp
    CB.snr_cut  = args.snr_cut
    CB.snr_thr  = args.snr_thr
    CB.mps      = args.mps
    CB.duration = args.duration
    CB.classify_args = [args.snr_cut,args.mem_cut,args.dm_cut,args.nbeam_cut, \
                        args.filter_cut,args.beam_mask,args.chan1,args.chanbw, \
                        args.tsamp]
    CB.start_browser()
