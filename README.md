# twocover-descent

This is an implementation of two-cover descent for genus 2 curves with a rational Weierstrass point. The theory and formulas are explained in the following paper: [https://arxiv.org/abs/2009.10313](https://arxiv.org/abs/2009.10313)

The main files are `twocovers.m` and `twocover-processor.sage`, containing the [Magma](http://magma.maths.usyd.edu.au/magma/) implementation and a Sage script that handles data processing and control flow for the algorithm using Sage's interface to Magma. Here is an example of how the script is used:

```
sage twocover-processor.sage --database data/g2c_curves-r2-w1.json --label 6443.a.6443.1 --output_directory ./results --stages all
```

An incomplete [Sage](https://www.sagemath.org/) implementation can be found in `twocovers.sage`. Note that the Sage code only implements geometric constructions, since Sage doesn't have the same facilities for computing Mordell-Weil groups or running the elliptic Chabauty method.
