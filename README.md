# A nonlinear voice from GW250114 ringdown
Yi-Fan Wang <sup>1,2</sup>, Sizheng Ma <sup>3</sup>, Neev Khera <sup>4</sup>, Junquan Su<sup>4</sup>, Huan Yang <sup>4</sup>

<sub>1.Max-Planck-Institut für Gravitationsphysik (Albert-Einstein-Institut), Am Mühlenberg 1, D-14476 Potsdam, Germany</sub>   
<sub>2.Purple Mountain Observatory, Chinese Academy of Sciences, Nanjing 210034, China</sub>   
<sub>3.Perimeter Institute for Theoretical Physics, Waterloo, ON N2L2Y5, Canada</sub>  
<sub>4.Department of Astronomy, Tsinghua University, Beijing 100084, China</sub>  

## Introduction
Gravitational-wave astronomy, by detecting ripples in spacetime, has opened a new window to observe compact objects and probe theories of gravity in the nonlinear strong-field regime. The ringdown signal of a binary black hole merger contains a superposition of damped sinusoids known as quasi-normal modes [1], whose frequencies are completely determined by the mass and spin of the remnant black hole and form the basis of \textit{black hole spectroscopy} [2-4]. A crucial prediction yet to be observationally confirmed is the existence of quadratic quasi-normal modes, which represent fundamental properties associated with wave–wave coupling in general relativity, and the leading mode is predicted to be detectable with next-generation ground-based detectors [5-7] using traditional methods. Here we show the first observational evidence for a set of quadratic quasi-normal modes in the ringdown of the binary black hole merger GW250114, the loudest gravitational-wave event detected to date, enabled by a novel analysis. These nonlinear modes result from the quadratic coupling of the linear $(2,2,n)$ modes with $n\leq3$. Starting the analysis at a time corresponding to four times the remnant mass ($M_\mathrm{f}$) after the merger, the evidence for their presence reaches a Bayes factor of 62. A phenomenological test allowing these modes to deviate from the theoretical prediction rejects the zero-amplitude hypothesis at a significance of 3.4 $\sigma$, while the inferred amplitude and complex frequency are consistent with the prediction of general relativity. This finding provides the first observational evidence of gravitational wave-wave interaction and extends black hole spectroscopy from the linear to the nonlinear regime. It also establishes a new direction for testing the fundamental nonlinear structure of general relativity with the most extreme gravity.

## Results & Reproduction

We use the following softwares and data to perform this work: 
 - [`pycbc`](https://github.com/gwastro/pycbc): Core package to analyze gravitational-wave data, find signals, and study their parameters.
 - [`tgr`](https://github.com/yi-fan-wang/TestingGR_with_Gravwaves): a pycbc waveform plugin for nonlinear quadratic quasi-normal modes waveforms.

The folder structure in this repository is:
 - `config`: configuration files used by `pycbc_inference` to obtain posterior files
 - jupyter notebooks (such as `fig01_bayes_factor_snr.ipynb`) : reproduce all figures in the paper

To reproduce this work and run jupyter notebooks, follow the instructions below:

### requirements installation

```bash
python -m venv env
source env/bin/activate
pip install -r requirements.txt
```

### NRSur7dq4 waveform data

Generating NRSur7dq4 waveforms and evaluating the remnant fits requires two additional model-data files: `NRSur7dq4_v1.0.h5` and
`NRSur7dq4Remnant_v1.0.h5`. These files must be available to LALSimulation in addition to the Python packages installed through `requirements.txt`.

Run the following commands from the repository root:

```bash
mkdir -p data/lalsuite

curl -fL --retry 3 \
    -o data/lalsuite/NRSur7dq4_v1.0.h5 \
    https://dcc.ligo.org/public/0198/T2500012/004/NRSur7dq4_v1.0.h5

curl -fL --retry 3 \
    -o data/lalsuite/NRSur7dq4Remnant_v1.0.h5 \
    https://dcc.ligo.org/public/0198/T2500012/004/NRSur7dq4Remnant_v1.0.h5

export LAL_DATA_PATH="$PWD/data/lalsuite${LAL_DATA_PATH:+:$LAL_DATA_PATH}"
```

Set `LAL_DATA_PATH` in each new terminal session before launching Jupyter or PyCBC Inference. If these files are already installed elsewhere, use their containing directory instead.

### LaTeX requirements for plotting

The plotting notebooks use Matplotlib with `text.usetex=True` and Computer Modern fonts. A working LaTeX installation is required; `dvipng` is also needed for raster rendering in Jupyter. These system dependencies are not installed by `pip`.

For Ubuntu/Debian, an example installation command is:

```bash
sudo apt-get install texlive-latex-extra texlive-fonts-recommended cm-super dvipng
```

For macOS with Homebrew, a full TeX installation can be installed with:

```bash
brew install --cask mactex-no-gui
```

After installing MacTeX, open a new terminal. Verify that the tools are available before starting Jupyter:

```bash
latex --version
dvipng --version
```

For a lightweight preview without LaTeX, replace the corresponding settings in each notebook's plotting configuration with:

```python
"text.usetex": False,
"font.serif": ["DejaVu Serif"],
```

This changes the figure typography but does not change the numerical analysis.

### GWOSC strain data and local file paths

To rerun parameter estimation, download the H1 and L1 strain data for GW250114_082203 from the GWOSC event page [here](https://gwosc.org/eventapi/html/O4_Discovery_Papers/GW250114_082203/v1/). Select the 4096-second, 16 kHz GWF files, starting at GPS time 1420877824.

Place the downloaded files in `data/gwosc/` under the repository root. The filenames currently listed by GWOSC are:

- `H-H1_GWOSC_O4b_16KHZ_R1-1420877824-4096.gwf`
- `L-L1_GWOSC_O4b_16KHZ_R1-1420877824-4096.gwf`

The released configuration files contain absolute paths from the original computing environment. Before running inference, update `frame-files` in BOTH `[inspiral__data]` and `[ringdown__data]` in each configuration you intend to use.

For the filenames above, use the following entries in both sections, while retaining the other configuration settings:

```ini
frame-files = H1:data/gwosc/H-H1_GWOSC_O4b_16KHZ_R1-1420877824-4096.gwf L1:data/gwosc/L-L1_GWOSC_O4b_16KHZ_R1-1420877824-4096.gwf
channel-name = H1:GWOSC-16KHZ_R1_STRAIN L1:GWOSC-16KHZ_R1_STRAIN
```

Downloading the strain data is not necessary for notebooks that only read the bundled `.npz` plotting data.

### Quick demo: reproduce Figure 1

After installing the Python dependencies and setting up LaTeX as described above, run the following commands from the repository root:

```bash
source env/bin/activate
jupyter lab fig01_bayes_factor_snr.ipynb
```

In JupyterLab, select **Run → Run All Cells**.

The notebook reads the bundled dataset `data/fig01_bayes_factor_snr.npz` and reproduces the Bayes-factor and signal-to-noise-ratio plot. This demo does not require downloading GWOSC
strain or NRSur model data, or rerunning parameter estimation.

The figure is displayed in the notebook and saved to:

```text
figures/bayes_factor_snr_gaussian_bar1sig_arrow_drawmedian.pdf
```

Run other notebooks to reproduce other figures.

### run PyCBC Inference

The command below is used for parameter estimation and generating the posteriors, which is not needed to run the jupyter notebooks above. Anyway, the command line to launch a PyCBC Inference run is (this should use a Linux or Mac operation system): 
```
OMP_NUM_THREADS=1 \
pycbc_inference --verbose \
    --seed 123456 \
    --config-file config/inference-all10.ini \
    --output-file posterior/H1L1-INFERENCE_ALL10.hdf \
    --nprocesses 32 \
    --force
```

It takes several days to complete the runs using 32 CPU cores.

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
