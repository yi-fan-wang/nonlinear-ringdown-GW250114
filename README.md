# A nonlinear voice from GW250114 ringdown
Yi-Fan Wang <sup>1</sup>, Sizheng Ma <sup>2</sup>, Neev Khera <sup>3,4</sup>, Huan Yang <sup>4</sup>

<sub>1.Max-Planck-Institut für Gravitationsphysik (Albert-Einstein-Institut), Am Mühlenberg 1, D-14476 Potsdam, Germany</sub>  
<sub>2.Perimeter Institute for Theoretical Physics, Waterloo, ON N2L2Y5, Canada</sub>  
<sub>3.Department of Physics, University of Guelph, Guelph, Ontario, Canada N1G 2W1</sub>  
<sub>4.Department of Astronomy, Tsinghua University, Beijing 100084, China</sub>  

## Introduction
The detection of quadratic quasi-normal modes would provide a direct probe into black hole nonlinear perturbations. We report the first observational evidence of a set of quadratic quasi-normal modes in the gravitational-wave ringdown of a binary black hole merger. Analyzing the signal from GW250114, we detect six nonlinear modes from the quadratic coupling of the fundamental $(2,2,0)$ mode and its first two overtones. At 5 final mass ($M_\mathrm{f}$) after the merger, the evidence for these nonlinear modes reaches a Bayes factor of 74. To single out these contributions, we employ recent theoretical progress to compute the waveforms and subtract the corresponding nonlinear modes from a numerical relativity surrogate waveform. Our data analysis uses a novel method that incorporates inspiral-merger inference results as a highly constraining prior for the ringdown inference. We further perform a test allowing for phenomenological deviations for the theoretically predicted amplitudes of the quadratic modes. The results show that an amplitude of zero is excluded at $3.0~\sigma$ significance level, while the theoretical expectation is consistent with the inference. This detection marks a first step towards observationally characterizing nonlinear perturbations in the ringdown of a black hole.

## Paper

[Arxiv Preprint](https://arxiv.org/abs/2601.05734)

## Results & Reproduction
 - `posterior`: posterior samples
 - `config`: configuration files used by `pycbc_inference` to obtain posterior files

To reproduce this work, one needs to install 
 - [`pycbc`](https://github.com/gwastro/pycbc)(v2.10.0 or the main branch): Core package to analyze gravitational-wave data, find signals, and study their parameters.
 - [`pytgr`](https://github.com/yi-fan-wang/TestingGR_with_Gravwaves)(v1.0 or the main branch): a pycbc waveform plugin for nonlinear quadratic quasi-normal modes waveforms.
 - Download the GW250114 strain data from [GWOSC](https://gwosc.org/eventapi/html/O4_Discovery_Papers/GW250114_082203/v1/)

An example command line to launch a PyCBC Inference run (this should use a Linux or Mac operation system): 
```
OMP_NUM_THREADS=1 \
pycbc_inference --verbose \
    --seed 123456 \
    --config-file config/inference-all10.ini \
    --output-file posterior/H1L1-INFERENCE_ALL10.hdf \
    --nprocesses 32 \
    --force
```

It takes a few minutes to install the dependent softwares on a computer, and O(1) days to complete the runs using 32 CPU cores.

## License and Citation

![Creative Commons License](https://licensebuttons.net/l/by-sa/4.0/88x31.png "Creative Commons License")

This work is licensed under a [Creative Commons Attribution-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-sa/4.0/).

We encourage use of these data in derivative works. If you use the material provided here, please cite the paper using the reference:

```
@article{Wang:2026rev,
    author = "Wang, Yi-Fan and Ma, Sizheng and Khera, Neev and Yang, Huan",
    title = "{A nonlinear voice from GW250114 ringdown}",
    eprint = "2601.05734",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    reportNumber = "LIGO-P2500804",
    month = "1",
    year = "2026"
}
```
