1. Purpose
The folder article contains the file that creates every figure in the article "Selective advantage of redirected helping in a viscous population," conducted during my thesis. The file main.ipynb creates all the results in this project. This file computes the inclusive fitness change of redirected helping in a viscous population in an infinite-island/patch model.

2 .System Requirements
I use Julia language (Julia 1.10.3) on Windows 10, 64 bits, to run the computations and to plot the results. One library is needed: PlotlyJS.

3. Instructions to use:
To run and plot the results, run the julia file in a terminal: julia main.jl. In case you need to add the package PlotlyJS, use "import Pkg" and Pkg.add("PlotlyJS"). 

4. Authorship
Author: Alan Flatrès, PhD candidate, aflatres@uwo.ca, a project supervised by Pr. Geoff Wild, gwild@uwo.ca, at the University of Western Ontario.

5. License:
MIT License

Copyright (c) [2024] [Alan Flatrès]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
