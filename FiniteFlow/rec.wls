#!/usr/bin/env -S math -script

<< FiniteFlow`;

rec[file_] := Module[
   {expr, dummy, vars},
   Print["Reconstructing ", file];
   expr = Together[Get[file]];
   vars = Variables[expr];
   FFNewGraph["g", "in", vars];
   FFAlgRatFunEval["g", "rf", {"in"}, vars, {expr}];
   FFGraphOutput["g", "rf"];
   FFReconstructFunctionMod[
      "g",
      vars,
      "PrintDebugInfo" -> 1,
      "MaxDegree" -> 200
   ];
];

rec["../data/aajamp"];
rec["../data/coeff_prop_4l"];
