# A nonlinear voice from GW250114 ringdown
Yi-Fan Wang <sup>1,2</sup>, Sizheng Ma <sup>3</sup>, Neev Khera <sup>4</sup>, Junquan Su<sup>4</sup>, Huan Yang <sup>4</sup>

<sub>1.Max-Planck-Institut für Gravitationsphysik (Albert-Einstein-Institut), Am Mühlenberg 1, D-14476 Potsdam, Germany</sub>   
<sub>2.Purple Mountain Observatory, Chinese Academy of Sciences, Nanjing 210034, China</sub>   
<sub>3.Perimeter Institute for Theoretical Physics, Waterloo, ON N2L2Y5, Canada</sub>  
<sub>4.Department of Astronomy, Tsinghua University, Beijing 100084, China</sub>  

## Introduction
Gravitational-wave astronomy, by detecting ripples in spacetime, has opened a new window to observe compact objects and probe theories of gravity in the nonlinear strong-field regime. The ringdown signal of a binary black hole merger contains a superposition of damped sinusoids known as quasi-normal modes [1], whose frequencies are completely determined by the mass and spin of the remnant black hole and form the basis of \textit{black hole spectroscopy} [2-4]. A crucial prediction yet to be observationally confirmed is the existence of quadratic quasi-normal modes, which represent fundamental properties associated with wave–wave coupling in general relativity, and the leading mode is predicted to be detectable with next-generation ground-based detectors [5-7] using traditional methods. Here we show the first observational evidence for a set of quadratic quasi-normal modes in the ringdown of the binary black hole merger GW250114, the loudest gravitational-wave event detected to date, enabled by a novel analysis. These nonlinear modes result from the quadratic coupling of the linear $(2,2,n)$ modes with $n\leq3$. Starting the analysis at a time corresponding to four times the remnant mass ($M_\mathrm{f}$) after the merger, the evidence for their presence reaches a Bayes factor of 62. A phenomenological test allowing these modes to deviate from the theoretical prediction rejects the zero-amplitude hypothesis at a significance of 3.4 $\sigma$, while the inferred amplitude and complex frequency are consistent with the prediction of general relativity. This finding provides the first observational evidence of gravitational wave-wave interaction and extends black hole spectroscopy from the linear to the nonlinear regime. It also establishes a new direction for testing the fundamental nonlinear structure of general relativity with the most extreme gravity.

## Paper

[Arxiv Preprint](https://arxiv.org/abs/2601.05734)

## Results & Reproduction
 - `config`: configuration files used by `pycbc_inference` to obtain posterior files
 - jupyter notebooks: reproduce all figures in the paper
```bash
python -m venv env
source env/bin/activate
pip install -r requirements.txt
jupyter lab
```

We use the following softwares and data to perform this work: 
 - [`pycbc`](https://github.com/gwastro/pycbc)(v2.10.0 or the main branch): Core package to analyze gravitational-wave data, find signals, and study their parameters.
 - [`tgr`](https://github.com/yi-fan-wang/TestingGR_with_Gravwaves)(the main branch): a pycbc waveform plugin for nonlinear quadratic quasi-normal modes waveforms.
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

## Change Log

We have substantially revised the paper in a v2 of the arXiv submission.

## License and Citation

![Creative Commons License](https://licensebuttons.net/l/by-sa/4.0/88x31.png "Creative Commons License")

This work is licensed under a [Creative Commons Attribution-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-sa/4.0/).

We encourage use of these data in derivative works. If you use the material provided here, please cite the paper using the reference:

```
@article{Wang:2026rev,
    author = "Wang, Yi-Fan and Ma, Sizheng and Khera, Neev and Su, Junquan and Yang, Huan",
    title = "{A nonlinear voice from GW250114 ringdown}",
    eprint = "2601.05734",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    reportNumber = "LIGO-P2500804",
    month = "1",
    year = "2026"
}
```
