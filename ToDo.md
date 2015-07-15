BUGS:
- MemoryError: store events on the fly
- High-DM astropysical pulses are removed (es. L199854_SAP0_BEAM72)
- Sometime infinite loop in parallelization

PLOT:
- BUG: in beam plot: histograms don't include top_candidates (check also general plot)
- Cumulative SNR instead of counts/SNR
- Histogram bins smaller
- Histogram bins logarithmic or plot non logarithmic
- Highlight time spans affected by RFI in general plot
- Histogram also for only pulses brighter than threshold (or maybe cumulative SNR) 
- General plot: obs and sap names
- General plot: legend discrete 
- Pulses shape: name of sap and beam and best/top
- General and beam plots: vertical lines also on central hist
- Change dim/SNR scaling of points
- General plot, top-right: y-axis limits accordingly to actual data values
- Single pulses: plot DM and Time instead of DM_c and Time_c
- In beam plot: plot only SNR>6
- In beam plot: numbers in white over the stars

OUTPUT:
- Folder with general best plots
- Save events for each beam and close them

ANALYSIS:
- Better removal of affected time bins
- Going below DM5
- Substitute time in alignment and grouping with sample
- Beam comparison: time range ~0.05s
- Remove pulses at the end of time
- Parameter for time spans affected by RFI in the whole sap
- Study better time alignment with high-DM pulsar
- Modify filters for the aligned pulse shape
- Study other possible filters: fit shape, pulse simmetry, pulse straight
- Substitute duration filter with downfact filter
- Alert for pulsars (more than x pulses at one DM, more than y pulses brighter than S in one DM)
- Pulse spectra
- Maybe possible to consider only bright pulses that appear in more beams
- Tollerance on DM steps in grouping relative to SNR of the event
- Filters also in incoherent beams

GENERAL
- Write Readme
- Re-write Installation

FRB
- Multimoment analysis: write out dynamic spectrum
- Matched filtering: modify single_pulse_search.py