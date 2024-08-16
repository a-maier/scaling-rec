- `coeff_prop_4l` is taken from the differential equation for a
  four-loop propagator. The propagator has the edges `[[0,1,0],
  [2,0,0], [4,2,1], [4,5,0], [5,0,0], [1,6,1], [7,2,1], [3,4,1],
  [6,5,0], [6,7,1], [7,3,0]]`, where the first two entries denote the
  adjacent vertices and the third entry the mass. External legs are
  connected to vertices 1 and 3. The rational function in
  `coeff_prop_4l` is the coefficient of the basis integral with indices
  `1,0,1,0,1,0,1,1,0,0,0`.

- `aajamp` is taken from the results of
  [arXiv:2105.04585](https://arxiv.org/abs/2105.04585):
  In
  [R2_qg_CF2_LmmpE.m](https://gitlab.msu.edu/vmante/aajamp-symb/-/blob/master/helicity_remainders/2loop/R2_qg_CF2_LmppE.m)
  it is the largest argument of `ratio` by Mathematica's
  `ByteCount`. After selecting the coefficient, the values r1,...,r22
  are replaced by x23, x34, x45, x51 using
  [definitions.m](https://gitlab.msu.edu/vmante/aajamp-symb/-/blob/master/aux/definitions.m)
