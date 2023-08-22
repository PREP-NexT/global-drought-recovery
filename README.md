# global-drought-recovery
Python code to Quantify global drought recovery probability based on Vine Copula

## License
This work is licensed under the GNU General Public License v3.0. For details please read LICENSE.txt.

## Usage
The scripts should be run in the following order:

#### 1. Gloabl_drought_events.py
Loads the data and extract characteristics of drought events.

#### 2. Find_severe_droughts.py
Find the most severe drought event for each grid during the past 66 years (1951-2016)

#### 3. Global_drought_recovery_probability.py
Calculate the likelihoods of drought recovery in both historcial (1951-1983) and present (1984-2016) periods
Evaluate whether changes in recovery probability between historical and present periods are statistically significant

#### 4. Elasticity_analysis.py
Calculate the response of drought recovery probability to precipitation changes

#### 5. Plot
(1) Plot the global recovery probability
(2) Plot the relative changes in recovery probability between historical and present periods at subcontinent scales
(3) Plot the response of drought recovery probability to precipitation changes under various climate scenarios
