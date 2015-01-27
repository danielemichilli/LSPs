BUG:
- Sort top_candidates
- Top_candidates in inf and plot different
- Only 10 top_candidates in inc beams
- Rank in top_candidates.inf wronf without inc beam

- General plot: high-right plot, axis range must include top and best
- Single pulse plots: enlarge points dimension
- General plots: enlarge points dimension
- General plots: diminish square dimension
- Single pulse plots: report obs name, sigma max, rank number
- Histograms in general plot
- Histogram bins smaller
- Histogram bins logarithmic or plot non logarithmic
- Beam plots: in top-right plot move numbers to the left
- Beam plots: give numbers to squares
- Move outside grid in plots
- In top_candidates.inf incoherent pulses are reported twice
- Remove time spans affected by RFI in general plot

- Better study parameters of grouping
  - Better study filter values
- Substitute duration filter with downfact filter
- Divide the general plot for the three saps
- Plot horiontal lines in DM trial changes
- Plots of pulse shapes larger when few pulses

TO TEST:
- Plots smaller
- Single pulse plots: circle dimensions relative to max
- Histogram also for only pulses brighter than threshold
- Strongest pulse in other list
- Strongest pulse only for DM > x

- General plot with counts in every beam
- Alert for pulsars (more than x pulses at one DM, more than y pulses brighter than S in one DM)
- Pulse spectra
