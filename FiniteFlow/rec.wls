#!/usr/bin/env -S math -script

<< FiniteFlow`;

rec[file_] := Module[
   {expr, dummy},
   Print["Reconstructing ", file];
   expr = Get[file];
   FFSparseSolve[
      {dummy == Together[expr]},
      {dummy},
      "PrintDebugInfo" -> 1,
      "MaxDegree" -> 200,
      "MaxPrimes" -> 200
   ];
];

rec["../data/aajamp"];
rec["../data/coeff_prop_4l"];
