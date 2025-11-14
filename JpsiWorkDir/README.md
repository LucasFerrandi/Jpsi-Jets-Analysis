# JpsiWorkDir

This is where O2 task are run, such as `dqEfficiency` and `tableMaker`. Typical content is:

- txt files containing lists of input files
- `run<...>.sh` scripts containing a pipeline with a list of O2 tasks and needed helper tasks
- json files containing the configuration for those tasks
- `AnalysisResult<...>.root` containing histograms produced by tasks
- `AO2D<...>.root` containing ROOT trees produced by tasks (usually meant for further analysis)