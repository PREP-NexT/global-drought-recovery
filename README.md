# Probabilistic Assessment of Global Drought Recovery and Its Response to Precipitation Changes

This repository contains the necessary Python code to quantify global drought recovery probability based on Vine Copula, introduced in the [AGU 2024](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023GL106067) paper:

*<center>"Probabilistic Assessment of Global Drought Recovery and Its Response to Precipitation Changes".</center>*


</br>


<center>
<img alt="fig1" width="800px" src="fig1.png">
</center>

## Usage
The scripts should be run in the following order:

#### 1. `global_drought_events.py`
Loads the data and extract characteristics of drought events

#### 2. `find_severe_droughts.py`
Find the most severe drought event for each grid during the past 66 years (1951-2016)

#### 3. `global_drought_recovery_probability.py`
Calculate the likelihoods of drought recovery in both historcial (1951-1983) and present (1984-2016) periods
Evaluate whether changes in recovery probability between historical and present periods are statistically significant

#### 4. `elasticity_analysis.py`
Calculate the response of drought recovery probability to precipitation changes

#### 5. `plot.py`
(1) Plot the global recovery probability
(2) Plot the relative changes in recovery probability between historical and present periods at subcontinent scales
(3) Plot the response of drought recovery probability to precipitation changes under various climate scenarios

## Citing this work

```bibtex
@article{zhang2024probabilistic,
  title={Probabilistic assessment of global drought recovery and its response to precipitation changes},
  author={Zhang, Limin and Yuan, Fei and He, Xiaogang},
  journal={Geophysical Research Letters},
  volume={51},
  number={1},
  pages={e2023GL106067},
  year={2024},
  publisher={Wiley Online Library}
}
```

## Code License
This work is licensed under the GNU General Public License v3.0. For details please read LICENSE.txt.
