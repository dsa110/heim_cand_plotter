#!/usr/bin/env python

import os,sys
import numpy as np

class Cand_Classifier(object):
    def __init__(self, snr=6.5, mem=3, dm=1.5, nbeam=3, filt=10, \
                       beamm=(1<<7)-1, chan1=1510., chanbw=0.5859375, \
                       tsamp=54.613e-6):
        self.snr_cut     = snr
        self.members_cut = mem
        self.dm_cut      = dm
        self.nbeams_cut  = nbeam
        self.filter_cut  = filt
        self.beam_mask   = beamm
        self.chan_1      = chan1  #MHz
        self.chan_bw     = chanbw #MHz
        self.tsamp       = tsamp  #sec

    def is_too_dim(self, cand):
        return cand['snr'] < self.snr_cut

    def is_noise(self, cand):
        try:
            return cand['members'] < self.members_cut
        except ValueError:
            return False

    def is_lowdm_rfi(self, cand):
        return cand['dm'] < self.dm_cut

    def is_too_fat(self, cand):
        return cand['filter'] > self.filter_cut

    def is_too_narrow(self, cand):
        min_chan_smear = (8.3e3 * cand['dm'] * self.chan_bw)/(self.chan_1**3)
        width = np.array([2**x for x in cand['filter']],dtype='float')*self.tsamp
        return width < min_chan_smear

    def is_masked(self, beam):
        return ((1<<beam) & self.beam_mask) == 0

    def is_hidden(self, cand):
        return ( self.is_masked(cand['beam']) |
                 ((self.is_masked(cand['beam']) != True) &
                  (cand['beam'] != cand['prim_beam'])) )

    def is_coincident(self, cand):
        #Mark the primary candidates (with highest snr) as coinc_rfi
        snr_prim_coinc_cands = cand['max_snr'][cand['beam'] != cand['prim_beam']]
        snr_prim_coinc_cands = np.unique(snr_prim_coinc_cands)
        coinc_rfi = np.isin(cand['snr'],snr_prim_coinc_cands)
        #Make fancier by using beam_mask
        #Take only events that occured in adjacent beams
        return coinc_rfi | (cand['nbeams'] > self.nbeams_cut)

    def is_coinc_rfi(self, cand):
        return cand['beam_mask'] > 2**16

    def categories(self, cand, multibeam, verbose=False):
        # Filter candidates based on classifications
        if verbose:
            print("Classifying candidates...\n")
        categories = {}
        is_hidden = self.is_too_dim(cand)
        is_noise  = (is_hidden==False) & (self.is_noise(cand) | self.is_lowdm_rfi(cand))
        is_coinc  = np.full(len(is_hidden),False)

        if multibeam:
            is_hidden = is_hidden | self.is_hidden(cand)
            is_noise  = (is_hidden==False) & (is_noise  | self.is_coinc_rfi(cand))
            is_coinc  = (is_hidden==False) & (is_noise ==False) & self.is_coincident(cand)

        is_wrong_width = (is_hidden==False) & (is_noise ==False) & (is_coinc ==False) & (self.is_too_fat(cand) | self.is_too_narrow(cand))
        is_valid  = (is_hidden==False) & (is_noise ==False) & (is_coinc ==False) & (is_wrong_width ==False)

        categories["hidden"]      = cand[is_hidden]
        categories["noise"]       = cand[is_noise]
        categories["coinc"]       = cand[is_coinc]
        categories["wrong_width"] = cand[is_wrong_width]
        categories["valid"]       = cand[is_valid]

        if verbose:
            print("Classified %i as hidden\n" % len(categories["hidden"]))
            print("           %i as noise spikes\n" % len(categories["noise"]))
            print("           %i as coincident RFI\n" % len(categories["coinc"]))
            print("           %i as wrong width\n" % len(categories["wrong_width"]))
            print("           %i as valid candidates\n" % len(categories["valid"]))

        return categories


class TextOutput(object):
    def __init__(self, multibeam=False):
        self.multibeam = multibeam

    def print_html(self, data, names):
        for cat in names:
            print("\n%s candidates:" % cat)
            if len(data[cat]) > 0:
                sys.stdout.write("<table width='100%' border=1 cellpadding=4px cellspacing=4px>\n")
                if self.multibeam:
                    sys.stdout.write("<tr><th align=left>SNR</th><th align=left>Time</th><th align=left>DM</th><th align=left>Filter [ms]</th><th align=left>Beam</th></tr>\n")
                else:
                    sys.stdout.write("<tr><th align=left>SNR</th><th align=left>Time</th><th align=left>DM</th><th align=left>Filter [ms]</th></tr>\n")

                for (i, item) in enumerate(data[cat]['snr']):
                    if self.multibeam:
                        sys.stdout.write ("<tr>" + \
                                          "<td>" + str(data[cat]['snr'][i]) + "</td>" + \
                                          "<td>" + str(data[cat]['time'][i]) + "</td>" + \
                                          "<td>" + str(data[cat]['dm'][i]) + "</td>" + \
                                          "<td>" + str(0.064 * (2 **data[cat]['filter'][i])) + "</td>" + \
                                          "<td>" + str(data[cat]['prim_beam'][i]) + "</td>" + \
                                          "</tr>\n")
                    else:
                        sys.stdout.write ("<tr>" + \
                                          "<td>" + str(data[cat]['snr'][i]) + "</td>" + \
                                          "<td>" + str(data[cat]['time'][i]) + "</td>" + \
                                          "<td>" + str(data[cat]['dm'][i]) + "</td>" + \
                                          "<td>" + str(0.064 * (2 **data[cat]['filter'][i])) + "</td>" + \
                                          "</tr>\n")
                sys.stdout.write("</table>\n")

    def print_text(self, data, names):
        for cat in names:
            print("\n%s candidates:" % cat)
            if len(data[cat]) > 0:
                if self.multibeam:
                    sys.stdout.write ("%-10s %-10s %-10s %-10s %-10s %-10s\n" \
                                       % ('snr', 'time', 'samp_idx', 'dm', \
                                          'filter', 'prime_beam'))
                else:
                    sys.stdout.write ("%-10s %-10s %-10s %-10s %-10s\n" \
                                       % ('snr', 'time', 'samp_idx', 'dm', 'filter'))
                for (i, item) in enumerate(data[cat]['snr']):
                    if self.multibeam:
                        sys.stdout.write ("%-10s " % data[cat]['snr'][i] + \
                                          "%-10s " % data[cat]['time'][i] + \
                                          "%-10s " % data[cat]['samp_idx'][i] + \
                                          "%-10s " % data[cat]['dm'][i] + \
                                          "%-10s " % data[cat]['filter'][i] + \
                                          "%-10s " % data[cat]['prim_beam'][i] + "\n")
                    else:
                        sys.stdout.write ("%-10s " % data[cat]['snr'][i] + \
                                          "%-10s " % data[cat]['time'][i] + \
                                          "%-10s " % data[cat]['samp_idx'][i] + \
                                          "%-10s " % data[cat]['dm'][i] + \
                                          "%-10s " % data[cat]['filter'][i] + "\n")

    def print_xml(self, data, names):
        for cat in names:
            print("\n%s candidates:" % cat)
            # get indicie list for sorting via snr
            snr_sorted_indices = [i[0] for i in sorted(enumerate(data[cat]['snr']), key=lambda x:x[1],reverse=True)]

            cand_i = 0
            for i in snr_sorted_indices:
                cand_i += 1
                if self.multibeam:
                    sys.stdout.write ("<candidate snr='" + str(data[cat]['snr'][i]) + \
                                               "' time='" + str(data[cat]['time'][i]) + \
                                               "' dm='" + str(data[cat]['dm'][i]) + \
                                               "' samp_idx='" + str(data[cat]['samp_idx'][i]) + \
                                               "' filter='" + str(data[cat]['filter'][i]) + \
                                               "' prim_beam='" + str(data[cat]['prim_beam'][i] + 1) + "'/>\n")
                else:
                    sys.stdout.write ("<candidate snr='" + str(data[cat]['snr'][i]) + \
                                               "' time='" + str(data[cat]['time'][i]) + \
                                               "' dm='" + str(data[cat]['dm'][i]) + \
                                               "' samp_idx='" + str(data[cat]['samp_idx'][i]) + \
                                               "' filter='" + str(data[cat]['filter'][i]) + "'/>\n")


def load_candidates(filename, verbose=False):
    if os.path.getsize(filename) > 0:
        try:
            # multibeam test data
            print("--aa--")
            all_cands = \
                np.loadtxt(filename, ndmin=1,
                    dtype={'names': ('snr','samp_idx','time','filter',
                                     'dm_trial','dm','members','begin','end',
                                     'nbeams','beam_mask','prim_beam',
                                     'max_snr','beam'),
                           'formats': ('f4', 'i4', 'f4', 'i4',
                                       'i4', 'f4', 'i4', 'i4', 'i4',
                                       'i4', 'i4', 'i4',
                                       'f4', 'i4')})
        except IndexError:
            try:
                # single beam test data
                print("--a--")
                all_cands = \
                    np.loadtxt(filename, ndmin=1,
                        dtype={'names': ('snr','samp_idx','time','filter',
                                         'dm_trial','dm','members','begin','end'),
                               'formats': ('f4', 'i4', 'f4', 'i4',
                                           'i4', 'f4','i4','i4','i4')})
            except:
                try:
                    print("--b--")
                    # DSA data pre-coincidencing
                    # itime=samp_idx, mjds=time (different units), ibox=filter, ibeam=beam, idm=dm_trial
                    all_cands = \
                        np.loadtxt(filename, ndmin=1,
                            dtype={'names': ('snr','if','samp_idx','time',
                                             'filter','dm_trial','dm','beam'),
                                   'formats': ('f4', 'i4', 'i4', 'f4',
                                               'i4', 'i4', 'f4', 'i4')})
                except:
                    try:
                        # DSA data post-coincidencing
                        all_cands = \
                            np.loadtxt(filename, ndmin=1,
                                dtype={'names': ('snr','if','samp_idx','time',
                                                 'filter','dm_trial','dm','members','beam'),
                                       'formats': ('f4', 'i4', 'i4', 'f4',
                                                   'i4', 'i4','f4', 'i4',
                                                   'i4')})

                    except IndexError:
                        raise

                    else:
                        multibeam = False
                else:
                    multibeam = False
            else:
                multibeam = False
        else:
            multibeam = True

            # Adjust for 0-based indexing
            all_cands['prim_beam'] -= 1
            all_cands['beam'] -= 1
    else:
        if verbose:
            print("Candidate file is empty!")
        all_cands = np.empty(0,
                        dtype={'names': ('snr','samp_idx','time','filter',
                                         'dm_trial','dm','members','begin','end'),
                               'formats': ('f4', 'i4', 'f4', 'i4',
                                           'i4', 'f4', 'i4', 'i4', 'i4')})
        multibeam = False

    if verbose:
        if multibeam:
            print("Loaded %i candidates from multiple beams\n" % len(all_cands))
        else:
            print("Loaded %i candidates from one beam\n" % len(all_cands))

    return all_cands, multibeam


if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Extract valid candidates from .cand files.")
    parser.add_argument('-v','--verbose', action="store_true")
    parser.add_argument('-V','--very_verbose', action="store_true")
    parser.add_argument('--snr_cut', type=float, default=6.5)
    parser.add_argument('--members_cut', type=int, default=3)
    parser.add_argument('--dm_cut', type=float, default=1.5)
    parser.add_argument('--nbeams_cut', type=int, default=3)
    parser.add_argument('--filter_cut', type=int, default=10)
    parser.add_argument('--beam_mask', type=int, default=(1<<7)-1)
    parser.add_argument('--chan_1', type=float, default=1510.)
    parser.add_argument('--chan_bw', type=float, default=0.5859375)
    parser.add_argument('--tsamp', type=float, default=54.613e-6)

    parser.add_argument('--cand_list_xml', action="store_true")
    parser.add_argument('--cand_list_html', action="store_true")

    parser.add_argument('cands_file')

    args = parser.parse_args()

    filename = args.cands_file
    if not os.path.isfile(filename):
        raise IOError("%s does not exist!" % filename)

    if args.very_verbose:
        names = ["hidden","noise","coinc","wrong_width","valid"]
        verbose = True
    elif args.verbose:
        names = ["valid"]
        verbose = True
    else:
        names = ["valid"]
        verbose = False

    all_cands, multibeam = load_candidates(filename, verbose)

    CC = Cand_Classifier(args.snr_cut, args.members_cut, args.nbeams_cut, \
                         args.dm_cut, args.filter_cut, args.beam_mask, \
                         args.chan_1, args.chan_bw, args.tsamp)
    categories = CC.categories(all_cands, multibeam, verbose)

    text_output = TextOutput()

    if args.cand_list_html:
      text_output.print_html(categories,names)
    elif args.cand_list_xml:
      text_output.print_xml(categories,names)
    else:
      text_output.print_text(categories,names)

    if verbose:
      print("\nDone\n")
