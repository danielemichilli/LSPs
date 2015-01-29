BUGS:
- High-DM astropysical pulses are removed (es. L199854_SAP0_BEAM72)
- Top_candidates in inf are wrong
- Only 10 top_candidates in inc beams
- Rank in top_candidates.inf wrong without inc beam
- Ignore corrupted tar files
- Check dimension of squares, sometimes not related with sigma apparently


PLOT:
- General plot: high-right plot, axis range must include top and best
- General plots: enlarge points dimension and diminish square dimension
- Single pulse plots: enlarge points dimension
- Single pulse plots: report obs name, sigma max, rank number
- Histogram bins smaller
- Histogram bins logarithmic or plot non logarithmic
- Beam plots: in top-right plot move numbers to the left
- Beam plots: give numbers to squares
- Remove time spans affected by RFI in general plot
- Divide the general plot for the three saps
- Plots of pulse shapes larger when few pulses
- Plots smaller
- Single pulse plots: circle dimensions relative to max
- Histogram also for only pulses brighter than threshold
- General plot with counts in every beam


OUTPUT:
- Folder with general best plots
- Sometimes top_candidates reported in the output are not the best and the same reported in the plot (es. L232335_SAP0_BEAM24)
- SAP0_BEAM12 has repeated best_pulses (es. L232335)
- In top_candidates.inf incoherent pulses are reported twice


ANALYSIS:
- Add filter duration/sigma 
- Better study parameters of grouping
  - Better study filter values
- Substitute duration filter with downfact filter
- Strongest pulse in other list
- Strongest pulse only for DM > x
- Alert for pulsars (more than x pulses at one DM, more than y pulses brighter than S in one DM)
- Pulse spectra
