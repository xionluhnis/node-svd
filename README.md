node-svd
========

Node.js module for [SVDLIBC](http://tedlab.mit.edu/~dr/SVDLIBC/) providing fast *Singular-Value-Decomposition* with ease.

Changes
=======
  * **0.1.3** - Now works with node version 0.12 and above

  * **0.1.2** - changed the Init method to InitSVD because of a name conflict with the new version of node

How to install
==============

You need **node-gyp** to be available as a command within the environment. Then, assuming it is available, simply use:

```Bash
npm install node-svd
```

How to use
==========

The module simply provide a function svd of signature:
```Javascript
/**
 * @param A [row1, row2 ...] matrix to apply svd onto
 * @param dim the number of singular values to keep, 0 for all (default)
 * @param settings { U, V, debug }
 * @return {d, U, S, V}
 */
```

About the settings:
* **debug** is a number from 0, 1 (default) and 2 describing the verbosity
* **U** is a boolean setting whether to untranspose the result *Ut* to *U* (true, by default)
* **V** is a boolean setting whether to untranspose the result *Vt* to *V* (false, by default)

The untransposition is provided because SVDLIBC result provides *Ut* and *Vt*, which are the transposed *U* and *V*, thus of dimensions (d by m) and (d by n).
*Vt* may be wanted, but if you want to recompute the initial *A*, you need *U* and not *Ut*, thus the default values.
Note that it doesn't add any overhead.

Example
=======

The basic is really simple:
```Javascript
var svd = require('svd').svd;
var res = svd(myMatrix, dim, settings);
// use res.U, res.S and res.V
```

Test
====

The basic example given in 'test.js' computes the svd of the simple matrix:
```Javascript
 A = [
  [1, 2],
  [3, 4],
  [5, 6]
];
```

The resulting U, S and V are displayed, and then the matrix A is recomputed (brute-force implementation of matrix multiplication) by simply multiplying back the parts together.
Floating point errors put aside, the results are the same (at least for that case).

Test check with Matlab / Octave
===============================

Simply use :
```Matlab
A = [1, 2; 3, 4; 5, 6];
[u, s, v] = svd(A)
B = u * s * v
C = u(1:3, 1:2) * s(1:2, 1:2) * v
```

LICENSE
=======

The MIT License

Copyright 2012 Alexandre Kaspar xion.luhnis@gmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
