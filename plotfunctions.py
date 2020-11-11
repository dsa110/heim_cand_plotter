#!/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

def TimeDMPlot(ax,data,duration=None,snr_cut=6.5,snr_thr=25.,mps=30.,multibeam=False,axrange=True,axlabel=True):
    ax.cla() #Clear the current axes

    ### Set global plot window parameters
    max_point_size = mps #diameter in points

    if axrange:
        DMS = np.array([])
        TT  = np.array([])
        for val in data.values():
            DMS = np.append(DMS,val['dm'])
            TT  = np.append(TT,val['time'])

        if duration != None:
            ax.set_xlim(0.,duration) #duration in seconds
        else:
            ax.set_xlim(0.,TT.max()*1.05)
        ax.set_ylim(1.,DMS.max()*1.05)
        ax.set_yscale('log')

    if axlabel == True:
        ax.set_xlabel('$\\rm Time\; (sec)$', size=12)
        ax.set_ylabel('$\\rm DM\;(pc\;cm^{-3})$', size=12)
    else:
        ax.tick_params(axis='x', labelbottom=False)

    # Set ticks and grid
    ax.grid(axis='y', which='both', linewidth=0.2)
    ax.set_axisbelow(True)
    yformat = plt.FormatStrFormatter('%i')
    ax.yaxis.set_major_formatter(yformat)

    # Plot data
    if (len(data['noise']['snr']) > 0):
        clipped_snr = np.clip(data['noise']['snr'],0.,snr_thr)
        norm_snr    = (clipped_snr-snr_cut)/(snr_thr-snr_cut)
        ax.scatter(data['noise']['time'],data['noise']['dm'], \
                   s=(max_point_size*norm_snr)**2, \
                   c='dimgray', alpha=0.5, marker='x', picker=5)

    if (len(data['coinc']['snr']) > 0):
        clipped_snr = np.clip(data['coinc']['snr'],0.,snr_thr)
        norm_snr    = (clipped_snr-snr_cut)/(snr_thr-snr_cut)
        ax.scatter(data['coinc']['time'],data['coinc']['dm'], \
                   s=(max_point_size*norm_snr)**2, \
                   c='slategrey', alpha=0.5, marker='+', picker=5)

    if (len(data['wrong_width']['snr']) > 0):
        clipped_snr = np.clip(data['wrong_width']['snr'],0.,snr_thr)
        norm_snr    = (clipped_snr-snr_cut)/(snr_thr-snr_cut)
        ax.scatter(data['wrong_width']['time'],data['wrong_width']['dm'], \
                   s=(max_point_size*norm_snr)**2, \
                   c=data['wrong_width']['filter'], cmap='viridis', \
                   vmin=0., vmax=12, alpha=0.25, marker='s', picker=5)

        if multibeam:
            for idx,time in enumerate(data['wrong_width']['time']):
                ax.text(time,data['wrong_width']['dm'][idx], \
                        data['wrong_width']['prim_beam'][idx].astype('str'), \
                        horizontalalignment='center', verticalalignment='center', \
                        size=6, color='dimgray', alpha=0.5)

    if (len(data['valid']['snr']) > 0):
        clipped_snr = np.clip(data['valid']['snr'],0.,snr_thr)
        norm_snr    = (clipped_snr-snr_cut)/(snr_thr-snr_cut)
        ax.scatter(data['valid']['time'],data['valid']['dm'], \
                   s=(max_point_size*norm_snr)**2, \
                   c=data['valid']['filter'], cmap='viridis', \
                   vmin=0., vmax=12, marker='o', picker=5)

        if multibeam:
            for idx,time in enumerate(data['valid']['time']):
                ax.text(time,data['valid']['dm'][idx], \
                        data['valid']['prim_beam'][idx].astype('str'), size=6, \
                        horizontalalignment='center', verticalalignment='center')



def DMSNRPlot(ax,ax_cbar,data,tsamp=54.613e-6):
    ax.cla() #Clear the current axes
    ax_cbar.cla()

    ### Set global plot window parameters
    max_snr = 100.

    DMS  = np.array([])
    SNRS = np.array([])
    for val in data.values():
        DMS  = np.append(DMS,val['dm'])
        SNRS = np.append(SNRS,val['snr'])

    ax.set_xlim(1.,DMS.max()*1.05)
    ax.set_xscale('log')
    ax.set_ylim(6.,max_snr)
    ax.set_yscale('log')

    # Set ticks and grid
    ax.grid(axis='x', which='both', linewidth=0.2)
    ax.set_axisbelow(True)
    xformat = plt.FormatStrFormatter('%i')
    ax.xaxis.set_major_formatter(xformat)

    yticks     = [6.,8.,10.,20.,40.,100.]
    yticklabel = ['%i' % tick for tick in yticks]
    ax.yaxis.set_ticks(yticks)
    ax.yaxis.set_ticklabels(yticklabel)
    ax.yaxis.set_minor_locator(plt.NullLocator())

    if SNRS.max() > max_snr:
        print("\n    WARNING: The SNR of a candidate is higher than the maximum SNRs plotted!\n")

    ax.set_xlabel('$\\rm DM\;(pc\;cm^{-3})$', size=12)
    ax.set_ylabel('$\\rm SNR$', size=12)

    # Plot data
    if (len(data['noise']['snr']) > 0):
        ax.scatter(data['noise']['dm'],np.clip(data['noise']['snr'],0.,max_snr), \
                   c='dimgray', alpha=0.5, marker='x', picker=5)

    if (len(data['coinc']['snr']) > 0):
        ax.scatter(data['coinc']['dm'],np.clip(data['coinc']['snr'],0.,max_snr), \
                   c='slategrey', alpha=0.5, marker='+', picker=5)

    if (len(data['wrong_width']['snr']) > 0):
        colormap = ax.scatter(data['wrong_width']['dm'],np.clip(data['wrong_width']['snr'],0.,max_snr), \
            c=data['wrong_width']['filter'], cmap='viridis', \
            vmin=0., vmax=12, alpha=0.25, marker='s', picker=5)

    if (len(data['valid']['snr']) > 0):
        colormap = ax.scatter(data['valid']['dm'],np.clip(data['valid']['snr'],0.,max_snr), \
            c=data['valid']['filter'], cmap='viridis', \
            vmin=0., vmax=12, marker='o', picker=5)

    # Add colorbar
    try: colormap
    except NameError:
        ax_cbar.axis("off")
        pass
    else:
        ax_cbar.axis("on")
        cbar = plt.colorbar(colormap,cax=ax_cbar,use_gridspec=True,label='$\\rm Boxcar width\;(ms)$')
        cticks = np.array(cbar.get_ticks())
        cbar.ax.set_yticklabels([x[0:5] for x in (2.**cticks * tsamp*1000.).astype('str')])
        cbar.set_alpha(1)
        cbar.draw_all()


def DMHistPlot(ax,all_cands,nbins=30,multibeam=False):
    #NOTE: this function takes in a list of all candidates,
    #      not the categorised list!
    ax.cla() #Clear the current axes

    ### Plot data
    dm_min = 1.
    if multibeam:
        for beam in range(all_cands['beam'].max()+1): #due to zero based indexing
            if beam in all_cands['beam']:
                cands = all_cands[all_cands['beam'] == beam]
                logbins = np.logspace(np.log10(dm_min),np.log10(cands['dm'].max()),nbins+1)
                vals,edges = np.histogram(cands['dm'],bins=logbins)
                ax.step(edges,np.append(vals,0.),where='post',label=str(beam))
                ax.legend(loc=9,ncol=4,fontsize=8)
    else:
        logbins = np.logspace(np.log10(dm_min),np.log10(all_cands['dm'].max()),nbins+1)
        ax.hist(all_cands['dm'],bins=logbins,histtype='step')

    ### Set global plot window parameters
    ax.set_xlim(dm_min,all_cands['dm'].max())
    ax.set_xscale('log')
    ax.set_ylim(1.,2000.) #As in heimdall plotter
    ax.set_yscale('log')

    # Set ticks and grid
    ax.grid(axis='x', which='both', linewidth=0.2)
    ax.set_axisbelow(True)
    formatter = plt.FormatStrFormatter('%i')
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)

    ax.set_xlabel('$\\rm DM\;(pc\;cm^{-3})$', size=12)
    ax.set_ylabel('$\\rm Candidate count$', size=12)


def DMHistPlot_multibeams(ax,all_cands,nbins=30):
    #NOTE: this function takes in a list of all candidates, not the categorised list!
    ax.cla() #Clear the current axes

    ### Plot data
    logbins = np.logspace(np.log10(1.),np.log10(all_cands['dm'].max()),nbins+1)
    ax.hist(all_cands['dm'],bins=logbins,histtype='step',orientation='horizontal')

    ### Set global plot window parameters
    ax.set_xlim(1.,2000.) #As in heimdall plotter
    ax.set_xscale('log')

    # Set ticks and grid
    ax.grid(axis='y', which='both', linewidth=0.2)
    ax.set_axisbelow(True)
    formatter = plt.FormatStrFormatter('%i')
    ax.xaxis.set_major_formatter(formatter)
    ax.tick_params(axis='y', labelleft=False)
    ax.tick_params(axis='x', bottom=False, labelbottom=False)


def DMSNRPlot_multibeams(ax,ax_cbar,data,tsamp=54.613e-6):
    ax.cla() #Clear the current axes
    ax_cbar.cla()

    ### Set global plot window parameters
    max_snr = 100.

    DMS  = np.array([])
    SNRS = np.array([])
    for val in data.values():
        DMS  = np.append(DMS,val['dm'])
        SNRS = np.append(SNRS,val['snr'])

    ax.set_xlim(6.,max_snr)
    ax.set_xscale('log')
    ax.set_xlabel('$\\rm SNR$', size=12)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.tick_params(axis='y', labelleft=False)

    # Set ticks and grid
    ax.grid(axis='y', which='both', linewidth=0.2)
    ax.set_axisbelow(True)
    xformat = plt.FormatStrFormatter('%i')
    ax.xaxis.set_major_formatter(xformat)

    xticks     = [6.,8.,10.,20.,40.,100.]
    xticklabel = ['%i' % tick for tick in xticks]
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_ticklabels(xticklabel)
    ax.xaxis.set_minor_locator(plt.NullLocator())

    if SNRS.max() > max_snr:
        print("\n    WARNING: The SNR of a candidate is higher than the maximum SNRs plotted!\n")

    # Plot data
    if (len(data['noise']['snr']) > 0):
        ax.scatter(np.clip(data['noise']['snr'],0.,max_snr),data['noise']['dm'], \
                   c='dimgray', alpha=0.5, marker='x', picker=5)

    if (len(data['coinc']['snr']) > 0):
        ax.scatter(np.clip(data['coinc']['snr'],0.,max_snr),data['coinc']['dm'], \
                   c='slategrey', alpha=0.5, marker='+', picker=5)

    if (len(data['wrong_width']['snr']) > 0):
        colormap = ax.scatter(np.clip(data['wrong_width']['snr'],0.,max_snr),data['wrong_width']['dm'], \
                   c=data['wrong_width']['filter'], cmap='viridis', \
                   vmin=0., vmax=12, alpha=0.25, marker='s', picker=5)

    if (len(data['valid']['snr']) > 0):
        colormap = ax.scatter(np.clip(data['valid']['snr'],0.,max_snr),data['valid']['dm'], \
                   c=data['valid']['filter'], cmap='viridis', \
                   vmin=0., vmax=12, marker='o', picker=5)

    # Add colorbar
    try: colormap
    except NameError:
        ax_cbar.axis("off")
        pass
    else:
        ax_cbar.axis("on")
        cbar = plt.colorbar(colormap,cax=ax_cbar,use_gridspec=True,label='$\\rm Boxcar width\;(ms)$')
        cticks = np.array(cbar.get_ticks())
        cbar.ax.set_yticklabels([x[0:5] for x in (2.**cticks * tsamp*1000.).astype('str')])
        cbar.set_alpha(1)
        cbar.draw_all()

def TimeBeamPlot(ax,data,duration=None,snr_cut=6.5,snr_thr=25.,mps=30.,multibeam=False,axrange=True,axlabel=True):
    ax.cla() #Clear the current axes

    ### Set global plot window parameters
    max_point_size = mps #diameter in points

    if axrange:
        B = np.array([])
        TT  = np.array([])
        for val in data.values():
            B = np.append(B,val['beam'])
            TT  = np.append(TT,val['time'])

        if duration != None:
            ax.set_xlim(0.,duration) #duration in seconds
        else:
            ax.set_xlim(0.,TT.max()*1.05)
        ax.set_ylim(1.,B.max()*1.05)

    if axlabel == True:
        ax.set_xlabel('$\\rm Time\; (sec)$', size=12)
        ax.set_ylabel('$\\rm Beam\;(pc\;cm^{-3})$', size=12)
    else:
        ax.tick_params(axis='x', labelbottom=False)

    # Set ticks and grid
    ax.grid(axis='y', which='both', linewidth=0.2)
    ax.set_axisbelow(True)
    yformat = plt.FormatStrFormatter('%i')
    ax.yaxis.set_major_formatter(yformat)

    # Plot data
    if (len(data['noise']['snr']) > 0):
        clipped_snr = np.clip(data['noise']['snr'],0.,snr_thr)
        norm_snr    = (clipped_snr-snr_cut)/(snr_thr-snr_cut)
        ax.scatter(data['noise']['time'],data['noise']['beam'], \
                   s=(max_point_size*norm_snr)**2, \
                   c='dimgray', alpha=0.5, marker='x', picker=5)

    if (len(data['coinc']['snr']) > 0):
        clipped_snr = np.clip(data['coinc']['snr'],0.,snr_thr)
        norm_snr    = (clipped_snr-snr_cut)/(snr_thr-snr_cut)
        ax.scatter(data['coinc']['time'],data['coinc']['beam'], \
                   s=(max_point_size*norm_snr)**2, \
                   c='slategrey', alpha=0.5, marker='+', picker=5)

    if (len(data['wrong_width']['snr']) > 0):
        clipped_snr = np.clip(data['wrong_width']['snr'],0.,snr_thr)
        norm_snr    = (clipped_snr-snr_cut)/(snr_thr-snr_cut)
        ax.scatter(data['wrong_width']['time'],data['wrong_width']['beam'], \
                   s=(max_point_size*norm_snr)**2, \
                   c=data['wrong_width']['filter'], cmap='viridis', \
                   vmin=0., vmax=12, alpha=0.25, marker='s', picker=5)

        if multibeam:
            for idx,time in enumerate(data['wrong_width']['time']):
                ax.text(time,data['wrong_width']['beam'][idx], \
                        data['wrong_width']['prim_beam'][idx].astype('str'), \
                        horizontalalignment='center', verticalalignment='center', \
                        size=6, color='dimgray', alpha=0.5)

    if (len(data['valid']['snr']) > 0):
        clipped_snr = np.clip(data['valid']['snr'],0.,snr_thr)
        norm_snr    = (clipped_snr-snr_cut)/(snr_thr-snr_cut)
        ax.scatter(data['valid']['time'],data['valid']['beam'], \
                   s=(max_point_size*norm_snr)**2, \
                   c=data['valid']['filter'], cmap='viridis', \
                   vmin=0., vmax=12, marker='o', picker=5)

        if multibeam:
            for idx,time in enumerate(data['valid']['time']):
                ax.text(time,data['valid']['beam'][idx], \
                        data['valid']['prim_beam'][idx].astype('str'), size=6, \
                        horizontalalignment='center', verticalalignment='center')
