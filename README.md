# A Fiery Death and State-Based Estimation
This project deals with the so-called Kalman filter, an algorithm for a recursive solution to the discrete- data linear filtering problem. Originally published in a 1960 paper by R.E. Kalman, the algorithm provides an efficient recursive computational means to estimate the state of a process, in a way that minimises the mean of the squared error (or variance).

## Build and run
Make sure you have cpgplot installed, the C-callable version of pgplot <http://www.astro.caltech.edu/~tjp/pgplot/>

### Command line
Build:
```
make proj_5
```

Then run:
```
./proj_5
```

### Within CLion, JetBrains' IDE for C projects
The `CMakeLists.txt` file should give you everything you need to specify build targets and required libraries, allowing you to build and run within the IDE.

## Data files
You will need to adjust the path to the data source files `x_data_p06.dat` and `y_data_p06.dat` to provide a full absolute path to these from the root of your local filesystem.

## Preparation of report
You can build the finished pdf of the report using the `Project5.tex` file and your favourite LaTeX compilation program. I used TexStudio on Mac <http://www.texstudio.org>