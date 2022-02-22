# twocover-descent

This is an implementation of two-cover descent for genus 2 curves with a rational Weierstrass point. The theory and formulas are explained in the following paper: [https://arxiv.org/abs/2009.10313](https://arxiv.org/abs/2009.10313). The raw data for the results of running this code on a dataset of 7692 curves of Mordell–Weil rank 2 or 3 can be found here: [https://github.com/HastD/twocover-results](https://github.com/HastD/twocover-results)

The main files are `twocovers.m` and `twocover-processor.sage`, containing the [Magma](http://magma.maths.usyd.edu.au/magma/) implementation and a [SageMath](https://www.sagemath.org/) script that handles data processing and control flow for the algorithm using Sage's interface to Magma. Here is an example of how the script is used:

```
sage twocover-processor.sage --label 6443.a.6443.1 --output_directory ./results --stages all
```

An incomplete pure Sage implementation can be found in `twocovers.sage`. Note that the Sage code only implements geometric constructions, since Sage doesn't have the same facilities for computing Mordell-Weil groups or running the elliptic Chabauty method.

The code that generates the data for the "Examples" section of the paper is in `paper-examples.sh`. If you're running this on your own system, delete the lines starting with `module avail` (these are instructions for the BU Shared Computing Cluster, which I used to run the code) and set `$outdir` to the directory where you want the output files to be written. You'll also first need to run
```
sage --preparse twocover-processor.sage
```
to generate the `.sage.py` file; it's more efficient (and avoids possible I/O errors if running several processes concurrently) to preparse the Sage file once at the outset and then call Sage with the `--python` option on the generated Python file afterward.

The `--index` option can be used in place of `--label` to specify the label at the given index in the supplied list of labels (default `data/labels.json`; can be set to a different file using the `--label_list` option). This is convenient for running the code on a large list of curves.

The file `data/g2c_curves-r2-w1.json` contains the coefficients of the hyperelliptic polynomials of each curve in the aforementioned dataset of 7692 genus 2 curves. If you want to run the code on curves not in this dataset, you will need to create a new JSON file with the coefficients of the curves in an associative array in the following format:
```
{"label1": [[coefficients of f], [coefficients of h]], "label2": ...}
```
where the labels will be the inputs to the `--label` option and the curve has equation `y^2 + h(x)*y = f(x)`. Use the `--database` option to tell the script where to find the coefficients.
