#!/usr/bin/env python3
#
# Taken from:
# https://github.com/devanshkv/pysigproc
# Changed by Leon Houben to make it more minimalistic for use
# with "show_candidates.py"

import numpy as np
from pysigproc import SigprocFile


def pad_along_axis(array, target_length, loc='end', axis=0, \
                   pad_mode='median', **kwargs):
    """

    :param array: Input array to pad
    :param target_length: Required length of the axis or tuple with pad sizes (before,end)
    :param loc: Location to pad: start: pad in beginning, end: pad in end, else: pad equally on both sides
    :param axis: Axis to pad along
    :return:
    """
    axis_nb = len(array.shape)
    npad = [(0, 0) for x in range(axis_nb)]

    if type(target_length) == tuple:
        npad[axis] = target_length
    else:
        pad_size = target_length - array.shape[axis]
        if pad_size < 0:
            return array

        if loc == 'start':
            npad[axis] = (pad_size, 0)
        elif loc == 'end':
            npad[axis] = (0, pad_size)
        else:
            npad[axis] = (pad_size // 2, pad_size // 2)

    return np.pad(array, pad_width=npad, mode=pad_mode, **kwargs)


class Fil_Cand(SigprocFile):
    def __init__(self, fp=None, dm=None, tcand=None, MJDcand=None, \
                 fref=None, _filter=0, kill_mask=None, verbose=False):

        """
        :param fp: Filepath of the filterbank
        :param dm: DM of the candidate
        :param tcand: Time of the candidate in filterbank file (seconds)
        :param MJDcand: MJD of the candidate in filterbank file
        :param fref: Reference frequency used to calculate tcand (MHz)
        :param _filter: width of candidate in 2**filter
        :param min_samp: Minimum number of time samples to read
        :param kill_mask: Boolean mask of channels to kill
        """
        SigprocFile.__init__(self, fp)
        self.dm      = dm
        self.tcand   = tcand
        self.MJDcand = MJDcand
        self.fref    = fref
        self.width   = 2.**_filter #in number of samples
        self.data    = None
        self.kill_mask = kill_mask
        self.verbose   = verbose

        if self.tcand is None and self.MJDcand is None:
            raise ValueError("Either tcand or MJDcand must be given!")
        if self.tcand is None and self.MJDcand is not None:
            self.tcand = (self.MJDcand - self.tstart)*86400.
        if self.tcand is not None and self.MJDcand is None:
            self.MJDcand = self.tstart + self.tcand/86400.
        else:
            Tcand = (self.MJDcand - self.tstart)*86400.
            assert (abs(self.tcand-Tcand) <= self.tsamp), "Given tcand and MJDcand are NOT equal!"

        if self.fref is None:
            self.fref = self.freqs.max()

        # To reduce the influence of DM errors and needed data
        # set the reference freq to that freq from which the DM delay,
        # to the top and bottom of the band, is the same.
        #
        # Assuming top band to be highst freq
        assert (self.foff < 0.), "Top of band is not the highest frequency! Flip the band."
        for idx,freq in enumerate(self.freqs):
            delays = self.dispersion_delay(self.freqs,freq,self.dm)
            if abs(delays[0]) > abs(delays[-1]):
                sym_freq = self.freqs[idx-1]
                break

        if self.fref != sym_freq:
            self.tcand += self.dispersion_delay(sym_freq,self.fref,self.dm)
            self.fref = sym_freq
            if self.verbose:
                print("\nReference freq set to %f MHz" % self.fref)


    @property
    def times(self):
        if self.data is None:
            return 0.
        else:
            ncand = self.get_ncand()
            return self.data_tsamp*np.arange(-ncand,self.data.shape[1]-ncand)

    @property
    def freqs(self):
        if self.data is None:
            return self.fch1 + np.arange(self.nchans)*self.foff
        else:
            return self.fch1 + np.arange(self.data.shape[0])*self.data_foff

    def dispersion_delay(self, nu=None, nu_ref=None, dms=None):
        """
        Calculates the dispersion delay at a specified DM
        :param dms: DM value to get dispersion delay
        :return:
        """
        if nu is None:
            nu = np.min(self.freqs)
        if nu_ref is None:
            nu_ref = np.max(self.freqs)
        if dms is None:
            dms = self.dm
        if dms is None:
            return None
        else:
            return 4.148808e3 * dms * (1./(nu**2) - 1./(nu_ref**2))

    def rebin1D(self, a, factor):
        return a.reshape(a.shape[0]//factor, -1).mean(1)

    def rebin2D(self, chanf, binf, a=None):
        if a is None and self.data is not None:
            sh = self.data.shape[0]//chanf,chanf,self.data.shape[1]//binf,binf
            self.data        = self.data.reshape(sh).mean(-1).mean(1)
            self.data_tsamp *= binf
            self.data_foff  *= chanf
            return
        elif type(a) is tuple:
            # Assuming tuple to have the form (x,y,xy)
            assert (len(a) == 3), "Given tuple is of wrong form!"
            t  = self.rebin1D(a[0],binf)
            f  = self.rebin1D(a[1],chanf)
            sh = a[2].shape[0]//chanf,chanf,a[2].shape[1]//binf,binf
            a  = a[2].reshape(sh).mean(-1).mean(1)
            return (t,f,a)
        else:
            sh = a.shape[0]//chanf,chanf,a.shape[1]//binf,binf
            return a.reshape(sh).mean(-1).mean(1)

    def get_ncand(self, bins_cut=0):
        TT      = np.arange(self.data.shape[1])*self.data_tsamp
        TT_diff = np.absolute(TT-self.data_tcand)
        ncand   = np.where(TT_diff == TT_diff.min())[0][0] - bins_cut
        return ncand

    def get_chunk(self, tstart=None, tstop=None):
        """
        Read the data around the candidate from the filterbank
        :param tstart: Start time in the fiterbank, to read from
        :param tstop: End time in the filterbank, to read till
        :return:
        """
        if tstart is None:
            tstart = self.tcand                                                  \
                        - 0.25*self.dispersion_delay()                           \
                        - abs(self.dispersion_delay(self.freqs.max(),self.fref)) \
                        - self.width * self.tsamp
        if tstop is None:
            tstop = self.tcand                                                   \
                        + 0.25*self.dispersion_delay()                           \
                        + abs(self.dispersion_delay(self.freqs.min(),self.fref)) \
                        + self.width * self.tsamp

        nstart = int(tstart / self.tsamp)
        nstop  = int(tstop / self.tsamp)
        if (nstop - nstart) % 2 == 1:
            nstop += 1

        # Get a chunk size that is a power of 2 long
        nbins = nstop-nstart
        diffs = np.absolute(2**np.arange(30) - nbins)
        chunk_size = 2**np.where(diffs == diffs.min())[0][0]

        nstart = nstart - (chunk_size - nbins)//2
        nstop  = nstop  + (chunk_size - nbins)//2

        if self.verbose:
            print("Set chunk size: %d bins" % nbins)
            print("Chunk size in pow 2: %d bins" % chunk_size)
            print("Delay across band: %d bins" % (np.round(self.dispersion_delay()/self.tsamp).astype(int)))

        if (nstart < 0) and (nstop > self.nspec):
            self.data = self.get_data(nstart=0, nsamp=self.nspec)
            if self.verbose:
                print('WARNING: Median padding data as nstart < 0 and nstop > nspec')
            self.data = pad_along_axis(self.data, (abs(nstart),nstop-self.nspec), axis=1)
            self.pad=True
        elif (nstart < 0) and (nstop < self.nspec):
            self.data = self.get_data(nstart=0, nsamp=nstop)
            if self.verbose:
                print('WARNING: Median padding data as nstart < 0')
            self.data = pad_along_axis(self.data, nstop-nstart, loc='start', axis=1)
            self.pad=True
        elif (nstart > 0) and (nstop > self.nspec):
            self.data = self.get_data(nstart=nstart, nsamp=(self.nspec-nstart))
            if self.verbose:
                print('WARNING: Median padding data as nstop > nspec')
            self.data = pad_along_axis(self.data, nstop-nstart, loc='end', axis=1)
            self.pad=True
        else:
            self.data = self.get_data(nstart=nstart, nsamp=(nstop-nstart))
            self.pad=False

        # Set data chunk specific parameters
        self.data_MJD     = self.tstart + \
                          (tstart - ((chunk_size - nbins)//2)*self.tsamp)/86400.
        self.data_tcand   = self.tcand - (self.data_MJD - self.tstart)*86400.
        self.data_foff    = self.foff
        self.data_tsamp   = self.tsamp

        if self.data.shape[1] != (nstop-nstart):
            raise ValueError("Data chunk does not have the desired length!")
        if self.data.shape[1] < np.round(self.dispersion_delay()/self.tsamp).astype(int):
            raise ValueError("Data chunk does not contain the full dispersion sweep!")
            
        # Set masked channels to zero
        if self.kill_mask is not None and self.kill_mask.any():
            assert len(self.kill_mask) == self.data.shape[0]
            data_copy = self.data.copy()
            data_copy[self.kill_mask,:] = 0
            self.data = data_copy
            del data_copy

    def dedisperse(self, dms=None, cut=True):
        """
        Dedisperse Frequency time data at a specified DM
        :param dms: DM to dedisperse at
        :return: dedispersed, valid dynamic spectrum
        """
        if dms is None:
            dms = self.dm

        if self.data is not None:
            delay_times = self.dispersion_delay(self.freqs,self.fref,dms)
            delay_bins  = -1*np.round(delay_times / self.data_tsamp).astype('int64')
            if cut:
                sweep_bins  = np.round(self.dispersion_delay() / self.data_tsamp).astype('int64')
                total_bins  = self.data.shape[1]-sweep_bins
                if total_bins % 2 == 1:
                    total_bins += 1
                dedispersed = np.zeros((self.data.shape[0],total_bins), dtype=self.dtype)
            else:
                dedispersed = np.zeros(self.data.shape, dtype=self.dtype)

            for idx in range(self.data.shape[0]):
                if cut:
                    dedispersed[idx] = np.roll(self.data[idx],delay_bins[idx])[delay_bins[0]:delay_bins[0]+total_bins]
                else:
                    dedispersed[idx] = np.roll(self.data[idx],delay_bins[idx])
        else:
            dedispersed = None

        if cut:
            ncand = self.get_ncand(delay_bins[0])
            t = self.data_tsamp*np.arange(-ncand,dedispersed.shape[1]-ncand)
        else:
            t = self.times
        f = self.freqs

        # Return x,y,z values
        return (t,f,dedispersed)

    def dedisperse_ts(self, dms=None, cut=True):
        """
        Dedisperse Frequency time data at a specified DM and return a time series
        :param dms: DM to dedisperse at
        :param cut: Return only valid spectra not wrapped by dedispersion
        :return: time series
        """
        if dms is None:
            dms = self.dm

        if self.data is not None:
            delay_times = self.dispersion_delay(self.freqs,self.fref,dms)
            delay_bins  = -1*np.round(delay_times / self.data_tsamp).astype('int64')
            if cut:
                sweep_bins  = np.round(self.dispersion_delay()/self.data_tsamp).astype('int64')
                total_bins  = self.data.shape[1]-sweep_bins
                if total_bins % 2 == 1:
                    total_bins += 1
                ts = np.zeros(total_bins, dtype=np.float32)
            else:
                ts = np.zeros(self.data.shape[1], dtype=np.float32)

            for idx in range(self.data.shape[0]):
                if cut:
                    ts += np.roll(self.data[idx],delay_bins[idx])[delay_bins[0]:delay_bins[0]+total_bins]
                else:
                    ts += np.roll(self.data[idx],delay_bins[idx])
        else:
            ts = None

        if cut:
            ncand = self.get_ncand(delay_bins[0])
            t = self.data_tsamp*np.arange(-ncand,ts.shape[0]-ncand)
        else:
            t = self.times

        # Return x,y values
        return (t,ts)

    def dmtime(self, dmsteps=256):
        """
        Generates DM-time array of the candidate by dedispersing at adjacent DM values
        dmsteps: Number of DMs to dedisperse at
        :return: DM-time array
        """
        dm_list     = self.dm + np.linspace(-self.dm, self.dm, dmsteps)
        dm_reversed = np.flip(dm_list,axis=0)
        dmt = np.zeros((dmsteps, self.data.shape[1]), dtype=np.float32)
        for ii in range(len(dm_list)):
            dmt[ii] = self.dedisperse_ts(dms=dm_reversed[ii], cut=False)[1]

        t = self.times

        # Return x,y,z values
        return (t,dm_list,dmt)

    def save_2Darray(self,arr=None,name=None):
        if self.data is None:
            return

        if arr is None:
            t   = self.times
            f   = self.freqs
            arr = self.data
        elif type(arr) is tuple:
            t   = arr[0]
            f   = arr[1]
            arr = arr[2]
        else:
            raise ValueError("Can't recognise the data to be saved.")

        if name is None:
            name = self.fp.name[:-4]+"_MJD"+str(self.data_MJD)

        axis2save = np.concatenate((t,f))
        np.save(name+"_xy",axis2save)
        np.save(name+"_spectrum",arr)

    def get_snr(self, time_series=None):
        """
        Calculates the SNR of the candidate
        :param time_series: time series array to calculate the SNR of
        :return:
        """
        if time_series is None and self.dedispersed is None:
            return None
        if time_series is None:
            x = self.dedispersed.mean(1)
        else:
            x = time_series
        argmax = np.argmax(x)
        mask = np.ones(len(x), dtype=np.bool)
        mask[argmax - self.width // 2:argmax + self.width // 2] = 0
        x -= x[mask].mean()
        std = np.std(x[mask])
        return x.max() / std
