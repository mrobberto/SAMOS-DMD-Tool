# SAMOS-DMD Tools
 
22/12/2024 - Version 5.0:


22/12/2024 - Version 4.9:

1) ranking is now pushed out of the notebook, i.e. the table is exported to .csf file and has to be ranked using e.g. MSexcel.

2) rewritten the slit configuration:
- a rectangular trace is created for each target; the trace/target is accepted if it does not overlap with the previous ones on a 2d array (simulates actual image).
- the order of injected traces is generated multiple times, randomly scrambling the weight=1 sample; the configuration providing the maximum number of targets and weight 1 targets is selected as mask0.
- the process is repeated on the remaining targets, finding mask1, mask2 etc up to a limit (currently set to 5, i.e. mask4)

3) removed all cells commented out using %%script false --no-raise-error.