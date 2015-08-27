BUGS:

PLOT:
- Histogram bins smaller
- Histogram bins logarithmic or plot non logarithmic
- Highlight time spans affected by RFI in general plot
- General plot: obs and sap names
- Pulses shape: name of sap and beam

OUTPUT:

ANALYSIS:
- Multimoment analysis
- Better removal of affected time bins
- Substitute time in alignment and grouping with sample
- Beam comparison: time range ~0.05s
- Remove pulses at the end of time
- Parameter for time spans affected by RFI in the whole sap
- Tollerance on DM steps in grouping relative to SNR of the event
- Study time misalignment on fake raw dataset
- RFI filters: 
    1. Many pulses at the same time but different DM
    2. Sigma / Number of beams in which a pulse appears (consider brightest and weakest pulses)
    3. Pulses at the same time in distant beams (create matrix of the beam disposition)

GENERAL
- Write Readme
- Re-write Installation

FRB
- Multimoment analysis: write out dynamic spectrum
- Matched filtering: modify single_pulse_search.py



Dynamic spectrum
