BUGS:
- MemoryError: store events on the fly
- High-DM astropysical pulses are removed (es. L199854_SAP0_BEAM72)

PLOT:
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

OUTPUT:
- Folder with general best plots
- Same name to folders inside hdf5 to all observations

ANALYSIS:
- Remove pulses at the end of tim
- Parameter for time spans affected by RFI in the whole sap
- Study better time alignment with high-DM pulsar
- Modify filters for the aligned pulse shape
- Study other possible filters: fit shape, pulse simmetry, pulse straight
- Better study parameters of grouping
  - Better study filter values
- Substitute duration filter with downfact filter
- Alert for pulsars (more than x pulses at one DM, more than y pulses brighter than S in one DM)
- Pulse spectra
- Maybe possible to consider only bright pulses that appear in more beams

GENERAL
- Write Readme
- Re-write Installation