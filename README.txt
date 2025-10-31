README for "Low Cost Sensors Deployment for Lake Monitoring: A Joint Optimization Strategy Integrating Quality¨CQuantity Trade-off and Value of Information Theory"

Version: 1.0  
Last updated: 2025-10-31

1. Overview

This dataset and accompanying MATLAB code provide all materials required to reproduce the results presented in the research on hierarchical monitoring network optimization. The framework integrates Value of Information (VOI) theory with network centrality analysis to construct efficient multi-level monitoring networks.

The files include all calculation outputs, performance metrics, convergence analysis, and supporting data used in the comprehensive analysis. All computations can be reproduced using the provided `hierarchicalSelection.m` script.

2. File List and Description

File Name 	Description
hierarchicalSelection.m	Main MATLAB function implementing the hierarchical monitoring strategy algorithm with VOI-centrality integration
VOI-calculation-results.xlsx	VOI calculation results containing detailed Value of Information statistics including effective VOI counts, ranges (min-max), mean, median, standard deviation, and coefficient of variation across different error ranges (10%~20%, 20%~30%, 30%~40%) for multiple datasets (201408, 201501, 201505)
VOI-Centrality-results-for-each-level.xlsx	Hierarchical analysis results showing VOI and centrality distributions across different levels (Level1-Level6) with comprehensive statistical summaries for each hierarchical tier
Information flow efficiency results.xlsx	Inter-level information transfer metrics quantifying the efficiency of information flow between consecutive hierarchical levels (Level1¡úLevel2, Level2¡úLevel3, etc.) across different sample groups and error ranges
Performance evaluation.xlsx	Comparative performance assessment of different hierarchical configurations, including coverage percentages per level, cumulative coverage analysis, and comparisons with expert experience benchmarks
convergence analysis.xlsx	Data convergence tracking documenting the optimization process convergence across multiple grades (Grade 1-7) with iteration counts (50, 100, 200, 400, 500) for comprehensive convergence behavior analysis
MI.xlsx	Mutual Information ranges providing MI value distributions (min-max ranges) across different grades (1-7) and error scenarios, used for information theory-based validation
MI-CI.xlsx	Mutual Information-Confidence Intervals  showing the relationship between mutual information and centrality measures across different grades and error ranges

3. Software Requirements

- MATLAB R2021a or later
- No additional toolboxes required - uses base MATLAB functionality only
- Compatible with Windows, macOS, and Linux platforms

4. How to Reproduce the Results

1. Download all files including `.m` function and `.xlsx` data files
2. Place all files in the same MATLAB working directory
3. Run the main function in MATLAB command window:
   matlab
   [levels, centrality_dynamic, coverage_history] = hierarchicalSelection();