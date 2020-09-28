// Copyright (c) 2020 Kinh Nguyen

// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#include "utils.hpp"

struct kinputs {
  const double beta;
  const double durInf;
  const double durLat;
  const double durQua;
  const double durTrt;
  const double nInfected;
  const double epi_time;
  const double * rho;
  const double * p_death;
  const int    * impulse_v;
  boost2D_ptr  POP;
  boost2D_ptr  eff;
  boost3D_ptr  M;
  boost2D_ptr impulse_a;
  // 
  double gamma() {return 1 - exp(-1/durInf);};
  double alpha() {return 1 - exp(-1/durLat);};
  double quara() {return 1 - exp(-1/durQua);};
  double treat() {return 1 - exp(-1/durTrt);};
  // double quara() {return 1 - exp(-1/durQua);};
  kinputs(const SEXP& inputs) :
    beta       (*REAL(get_value(inputs, "beta"))),
    durInf     (*REAL(get_value(inputs, "durInf"))),
    durLat     (*REAL(get_value(inputs, "durLat"))),
    durQua     (*REAL(get_value(inputs, "durQua"))),
    durTrt     (*REAL(get_value(inputs, "durTrt"))),
    nInfected  (*REAL(get_value(inputs, "nInfected"))),
    epi_time   (*REAL(get_value(inputs, "epi_time"))),
    rho        ( REAL(get_value(inputs, "rho"))),
    p_death    ( REAL(get_value(inputs, "p_death"))),
    impulse_v  ( INTEGER(get_value(inputs, "impulse_v"))),
    POP        ( REAL(get_value(inputs, "POP")), get_dim_2D(inputs, "POP")),
    eff        ( REAL(get_value(inputs, "eff")), get_dim_2D(inputs, "eff")),
    M          ( REAL(get_value(inputs, "M")), get_dim_3D(inputs, "M")),
    impulse_a  ( REAL(get_value(inputs, "impulse_a")), get_dim_2D(inputs, "impulse_a"))
  {}
};

extern "C" SEXP ksimc(SEXP inputs) {
  int np = 0;
  kinputs ip(inputs);

  double gamma = ip.gamma();
  double alpha = ip.alpha();
  double quara = ip.quara();
  double treat = ip.treat();

  boost1D N_age = ip.POP[ indices[1][_all] ];
  boost1D p_age = ip.POP[ indices[2][_all] ];
  // double N = sum_vector(N_age);

  const int S = 0, E = 1, Is = 2, Ic = 3, R = 4, D = 5, Q=6, T=7, 
    lambda=8, N_agegr = 16, N_state = 9;
  int matsize = ip.epi_time * N_agegr * N_state;
  SEXP output = PROTECT(NEW_NUMERIC(matsize)); ++np;
  memset(REAL(output), 0, matsize * sizeof(double));
  boost3D_ptr o(REAL(output), extents[ip.epi_time][N_agegr][N_state]);
  for (int i = 0; i < N_agegr; ++i) {
    o[0][i][Ic] = ip.nInfected;
    o[0][i][S]  = N_age[i] - (o[0][i][E] + o[0][i][Ic] + o[0][i][Is] + o[0][i][R]);
  }
  double SE, EIc, EIs, IsR, IcT, IcQ, Q2R, T2D;
  for (int i = 0; i < ip.epi_time-1; ++i) { // day
    if (ip.impulse_v[i])
      for (int j = 0; j < N_agegr; ++j)
          o[i][j][E] = ip.impulse_a[i][j];
    for (int j = 0; j < N_agegr; ++j) {     // me
      for (int k = 0; k < N_agegr; ++k)     // other
        for (int l = 0; l < 4; ++l)         // home > work > school > other
          o[i][j][lambda] += 
            ip.beta * ip.M[l][k][j] * ip.eff[i][l] * ((o[i][k][Ic] + o[i][k][Is]) / N_age[k]);

      SE  = o[i][j][lambda]       * o[i][j][S];
      EIc = alpha *    ip.rho[j]  * o[i][j][E]; // more old went infectious
      EIs = alpha * (1-ip.rho[j]) * o[i][j][E]; // more young went asymp
      IcQ = gamma * o[i][j][Ic] * (1 - ip.p_death[j]); // more young when to Q then recovered
      IcT = gamma * o[i][j][Ic] *      ip.p_death[j]; // more old went to treat
      IsR = gamma * o[i][j][Is]; // not counted towards statistics
      T2D = treat * o[i][j][T]; // treat some day then died
      Q2R = quara * o[i][j][Q]; // spent time in quarantine
      
      o[i+1][j][S]  = o[i][j][S] - SE;
      o[i+1][j][E]  = o[i][j][E] + SE - EIc - EIs;
      o[i+1][j][Ic] = o[i][j][Ic]     + EIc - IcQ - IcT;
      o[i+1][j][Is] = o[i][j][Is]           + EIs - IsR;
      o[i+1][j][T]  = o[i][j][T]                  + IcT - T2D;
      o[i+1][j][Q]  = o[i][j][Q]            + IcQ - Q2R;
      o[i+1][j][R]  = o[i][j][R]                  + Q2R;
      o[i+1][j][D]  = o[i][j][D]                        + T2D;
    }
  }
  SEXP out_dim = PROTECT(NEW_INTEGER(3)); ++np;
  INTEGER(out_dim)[0] = N_state;
  INTEGER(out_dim)[1] = N_agegr;
  INTEGER(out_dim)[2] = ip.epi_time;
  SET_DIM(output, out_dim);
  UNPROTECT(np);
  return output;
}