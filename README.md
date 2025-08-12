# project-ieeg-info-encoding
**A rapid theta network mechanism for flexible information encoding**

Quickly adapting our behavior depends on mechanisms that open and close a metaphorical ‘gate’ to let only relevant information into working memory. Using direct brain recordings in neurosurgical patients (i.e., intracranial EEG/iEEG), we discovered neocortical brain networks that rapidly come in and out of sync to direct this gate – a role previously attributed to a subcortical region called the striatum. Findings reveal a neocortical gating mechanism and demonstrate how human behavior can emerge from interactions across brain networks. Press release by [Northwestern News](https://news.feinberg.northwestern.edu/2023/06/01/study-establishes-fluctuating-gating-mechanisms-supporting-flexible-behavior/).

Our research aims to uncover mechanistic explanations of the neural basis of human behavior, that is, move from where to how. Our goals are multifaceted: (1) advance fundamental science by discovering new knowledge using rigorous, reproducible methods; and (2) advance translational applications in neurotechnology, precision medicine, and product development that are grounded in rigorous science. In this project, we showcase how approaching neuroscience as data science - with neurophysiological, neuroimaging, and behavioral datasets representing diverse data types, and analytics spanning signal processing, real-time behavior prediction, and graph theory - can address outstanding questions about information processing. Understanding how the human brain processes information can inform the design of artificial agents that are not just capable, but responsible.

These scripts produce the results reported in:

Johnson, EL, Lin, JJ, King-Stephens, D, Weber, PB, Laxer, KD, Saez, I, Girgis, F, D’Esposito, M, Knight, RT, Badre, D. A rapid theta network mechanism for flexible information encoding. _Nature Communications_ 14 (2023). [DOI](https://doi.org/10.1038/s41467-023-38574-7)

Project data are available on [OSF](https://doi.org/10.17605/OSF.IO/RX2ZD). These data have been prepared for analysis using the [pipeline for iEEG data curation](https://github.com/elizljohnson-projects/pipeline-ieeg-data-curation.git).

Publications or other papers using these scripts and/or data should cite the publication above.

Software:
- MATLAB 9.7 (R2019b; last tested with R2024b)
- Fieldtrip (fieldtrip-20210301; last tested with fieldtrip-20250114) - [download](https://www.fieldtriptoolbox.org/download)

Notes:
- Open and run the batch script (run_batch) to run through all scripts and reproduce the published results.
- The subfunctions folder contains the following functions:
  - zbaseline_200: z-score per-trial time-series data on the baseline using statistical bootstrapping, randomly drawing 200 data points per permutation. This script is task-agnostic.
  - zcont_5: z-score per-channel continuous data on itself using statistical bootstrapping, randomly drawing 5 data points per permutation. This script is task-agnostic.
  - split_cf_cl: split the data by condition and trial type. This script is task-specific.
- Data distribution plots were created using raincloud plotting scripts (v2) for MATLAB - [download](https://github.com/RainCloudPlots/RainCloudPlots/tree/master/tutorial_matlab).
- Brain data plots were created using BrainNet Viewer (v1.7) for MATLAB - [download](https://www.nitrc.org/projects/bnv).
