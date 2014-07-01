#!/usr/bin/env python
'''
Find pulses, sort by SNR, look for double detections.

Note: built using the new single pulse grouping algorithm. This script will
in future be filtering out bad data on the fly, allowing it to be run on bad
messy data sets.
'''
import os
import sys
import optparse
import time

from ssps import candidate
from ssps import pulse
from ssps import pulsetrain
from ssps import diagnostic
from ssps.support import check_delays_option

# -----------------------------------------------------------------------------
# -- Settings for data reduction pipeline -------------------------------------

DMS_ADJACENT = 4



def ssps_grab(searchoutdir):
    t0 = time.time()
  
    odir = '.'
    # read data
    print '=' * 77
    print 'Processing %s' % searchoutdir
    print 'Looking for datafiles.'
    spr = candidate.SinglePulseReaderBase(searchoutdir, delays_file, 30,
                                          options.lo_dm, options.hi_dm)

    print 'DM range after reading: [%.2f, %.2f]' % (spr.dms[0], spr.dms[-1])

    if options.ndms is not None:
        ndms = len(spr.dms)
        if ndms < options.ndms:
            print 'Too few DM trials, aborting!'
            print 'Needed %d and found %d DM trials.' % (options.ndms, ndms)
            sys.exit(1)
    # determine the basename for the output files:
    if options.o is None:
        dm_str = '_DM%.2f' % spr.dms[0]
        basename = spr.md_map[spr.dms[0]].data_file[:-len(dm_str)]
    else:
        basename = options.o
    print 'Looking for pulses.'
    # Call the rewritten ssps candidate grouping algorithm to find single
    # pulses by combining candidates across DM trials.
    pulses = pulse.group(spr, DMS_ADJACENT)
    pulses = pulse.annotate_pulses(spr, pulses)

    # Find the bright pulses and the dim ones, sort them for SNR.
    # find pulsar/RRAT candidates
    print 'Sifting groups, single pulse detections, into bright and dim ones.'
    print '  Threshold (peak) snr for bright pulses %f .' % options.snr
    print '  Minimum number of candidates per group %d .' % \
        options.min_n_candidates_bright
    print '  Minimum no. of bright pulses to consider a DM interesting %d.' % \
        options.min_n_repeats
    pulse_trains, rejects = \
        pulsetrain.extract_pulsetrains(spr, pulses, options.snr,
                                       options.min_n_repeats,
                                       options.min_n_candidates_bright,
                                       options.n_dms_bright,
                                       options.tlo_dm, options.thi_dm)

    n_bright = sum(len(x) for x in pulse_trains)
    print '  Found %d bright pulses near %d trial DMs.' % \
        (n_bright, len(pulse_trains))
    print 'Adding dim pulses to the bright ones found already.'
    print '  Minimum number of candidates per group %d .' % \
        options.min_n_candidates_dim
    pulse_trains, rejects = \
        pulsetrain.add_dim_pulses(spr, pulse_trains, rejects,
                                  options.min_n_candidates_dim)
    n_dim = sum(len(x) for x in pulse_trains) - n_bright
    print '  Found %d dim pulses.' % n_dim

    for keepers in pulse_trains:
        print 'DM %.2f had %d detections with max(snr) = %.2f.' % \
            (keepers[0].dm, len(keepers), max(x.snr for x in keepers))

    print '''TOOK %.2f SECONDS.''' % (time.time() - t0)

    # Determine the output file names:
    basenames = []
    for i, keepers in enumerate(pulse_trains):
        fn = basename + ('_%04d' % i) + ('_DM%.2f' % (keepers[0].dm))
        if not options.overwrite:
            if os.path.exists(os.path.join(options.odir, fn + '.xml')) or \
                    os.path.exists(os.path.join(options.odir,
                                   fn + '.peaks.singlepulse')):
                print '\nERROR: Output files exist, refusing overwrite!'
                sys.exit(1)

        basenames.append(fn)

    # Plot everything:
    n_plots = len(pulse_trains)
    for i in range(n_plots):
        # Determine the full path to this diagnostic plot:
        fn = os.path.join(options.odir, basenames[i] + '.xml')
        print 'Writing plot to %s' % fn
        # Determine what the links in the diagnostic plot should point to:
        if options.uselinkplaceholder:
            PREVIOUS_LINK = 'PREVIOUS_PLACEHOLDER'
            NEXT_LINK = 'NEXT_PLACEHOLDER'
        else:  # Remember, we want relative links here!
            if i > 0:
                PREVIOUS_LINK = basenames[i - 1] + '.xml'
            else:
                PREVIOUS_LINK = ''
            if i < n_plots - 1:
                NEXT_LINK = basenames[i + 1] + '.xml'
            else:
                NEXT_LINK = ''
        # Make the plot:
        diagnostic.plot(fn, spr, options, pulse_trains, i,
                        rejects, NEXT_LINK, PREVIOUS_LINK, marker_dms)
        # Write a file with the peaks of each pulse in the pulsetrain (in the
        # same format as the PRESTO .singlepulse files).
        pulse.write_arrivals(os.path.join(options.odir,
                             basenames[i] + '.peaks.singlepulse'), pulse_trains[i])
    print 'Took %.2f seconds!' % (time.time() - t0)
    
    
    return
