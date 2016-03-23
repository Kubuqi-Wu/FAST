# Installation

### Installation requirements
- Matlab (FAST has been developped and tested on Matlab R2015a)
- One of the following solvers :
  * Linprog (from the optimization toolbox)
  * Gurobi (recommended)
  * Mosek
  * GLPK

### Installation instructions
1. Unzip this package, and put it anywhere you want
2. Add to the Matlab path :
  * The src folder of this toolbox
  * Gurobi or Mosek if you intend to use it. If not, linprog should do the job, at least for small size problems. GLPK is also supported, but we do not advice to use it.

Once this is done, check out examples in the examples folder.

Let us know if you experience difficulties with the installation.
In particular, if your version of Matlab is not recent enough and we use a function you don't have, let us know. We can try to remove this dependency.

Go [here](http://baemerick.be/fast/demo.php) for more information and for a quick start example.