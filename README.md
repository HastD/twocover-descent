# twocover-descent

This is an implementation of two-cover descent for genus 2 curves with a rational Weierstrass point. The theory and formulas are explained in the following paper: [https://arxiv.org/abs/2009.10313](https://arxiv.org/abs/2009.10313)

The main files are `twocovers.m` ([Magma](magma.maths.usyd.edu.au/) implementation) and `twocovers.sage` ([Sage](https://www.sagemath.org/) implementation). Note that the Sage code only implements the geometric constructions, since Sage doesn't have the same facilities for computing Mordell-Weil groups or running the elliptic Chabauty method.

`twocovers-procedure.m` shows an example of usage. `twocovers-test.m` has miscellaneous (somewhat disorganized) example code.

