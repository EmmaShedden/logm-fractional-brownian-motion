# Fractional Brownian motion

## Description
This project is a continuation of our work simulating standard Brownian motions, which you can find [here](https://gitlab.eecs.umich.edu/logm/wn23/fractional-brownian-motion/standard-brownian-motion).

Fractional Brownian motion (fBm) is a generalization of classical Brownian motions and was first proposed by Mandelbrot and Van Ness in 1968 [[1](https://epubs.siam.org/doi/10.1137/1010093)]. Whereas classical Brownian motions describe continuous stochastic processes with independence of increments, in general, fBm exhibits long-range independence between increments. This allows it to model a larger class of real-world phenomena, but also makes it much more difficult so simulate accurately and efficiently. Fractional Brownian motion has been applied in various fields such as [medical imaging and robotics](https://dlib.bc.edu/islandora/object/bc-ir:102098), we can even use fractional Brownian motion for [composing music](https://www.sciencedirect.com/science/article/pii/0167278989902200)!

This repository contains:
- Implementation of the [Wood-Chan simulation method](https://drive.google.com/file/d/1BEjP1AHJWwW1HtJDZcKPLzWJ1wXDoxcW/view) for fractional Brownian motions (Python).
- Implementation of the Fast Fourier Transform (FFT) and inverse FFT on powers of two, a useful algorithm which our simulation method relies on for efficiency over naive approaches (Python).
    - Randomized unit tests of this implementation against standard library implementations in `scikit` and `numpy`.
    - This is for learning purposes; the standard library functions are used in our simulation code.
- Implementation of the [Series Representation of fBm](http://math.univ-lille1.fr/~ayache/COURSE-WavFrac.pdf) (Python).

## Installation
See our installation instructions [here](https://gitlab.eecs.umich.edu/logm/wn23/fractional-brownian-motion/standard-brownian-motion#installation).

## Usage
Usage instructions for all programs within this repository. Make sure you have completed the relevant installation for the files you want to run.

All instructions begin in the command line at the top level of the repo.

Depending on your setup, you may need to replace the command `python` with `python3` when running Python scripts.

### Fast Fourier Transform (FFT)
1. Run `cd fft/`.
2. Run `python main.py`. This runs randomized unit tests of our FFT implementation. The output to the terminal describes the test suite.
3. To modify the unit test configuration, open `main.py` in your editor and scroll down to the definition of `main()`.
    - Change `m` to run a different number of tests per input length. This must be the same for all input lengths.
    - Change `p` to test a different range of array lengths. All inputs have lengths equal to a power of two, from $1$ to $2^{p-1}$ inclusive.
    - To add specific hardcoded unit tests, pass them into `test()` per the example in the code.
    - For larger input arrays, it is necessary to increase the tolerance for machine error in floating point comparisons (`tol`).

### Wood-Chan Algorithm
1. Run `cd wood-chan/`
2. Run `python main.py`. This will generate plots for paths of fBm on the unit interval, discretized to $2^10 + 1$ timestamps. The $H$ values are approximately $\frac{1}{4}$, $\frac{1}{2}$, $\frac{3}{4}$, and $1$. The output is written to `wood-chan/fBm.png`.
3. To simulate different $H$ values or change the granularity of the time interval, open `main.py` in your editor and scroll down to the definition of `main()`.
    - Change `epsilon` to sample different $H$ values.
    - Change `q` to change the granularity of discretization of the time interval.

#### Example output
<p align="center"><img width="600" src="/wood-chan/fBm_9.png"></p>

### Series Representation of fBm
1. Run `cd series-rep/`
2. Refer to the instructions for your choice of installation in order to open `main.ipynb`.
3. Change any settings you wish in the cell at the top of the file. You can choose whether to display and/or store the output graphs, and whether to run intermediate tests of individual functions.
4. In the top menu selection, find __Kernel__ and select __Restart & Run All__.
    - Note that the generation of the images takes several minutes. This is normal.
5. If you choose to store the output graphs, they can be found in `[path-to-repo]/series-rep/`. You can view them in your IDE or any image viewer.
    - Extra information on how to interpret output graphs can be found in inline comments in the code.

#### Example output
<p align="center"><img width="600" src="/series-rep/fBm_4.png"></p>

## Support
If you encounter trouble running our code, you can email the following:
- xiaoranc@umich.edu
- lesteven@umich.edu
- emshedde@umich.edu

## Authors and acknowledgment
__Authors:__ Xiaoran Chen, Seok Jin Lee, & Emma Shedden

__Mentor:__ Hai Le

__Instructors:__ Ahmad Barhoumi, Sam Hansen

__Institution:__ University of Michigan

## Project status
Complete.
