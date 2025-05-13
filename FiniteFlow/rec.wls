#!/usr/bin/env -S math -script

<< FiniteFlow`;

rec[file_, args___] := Module[
   {expr, dummy, vars},
   expr = Together[Get[file]];
   vars = Variables[expr];
   FFNewGraph["g", "in", vars];
   FFAlgRatFunEval["g", "rf", {"in"}, vars, {expr}];
   FFGraphOutput["g", "rf"];
   FFReconstructFunctionMod[
      "g",
      vars,
      "PrintDebugInfo" -> 1,
      "MaxDegree" -> 200,
      args
   ];
];

Do[
   If[Length[opt] > 0,
      Print["Reconstructing ", file, " with options ", opt],
      Print["Reconstructing ", file]
   ];
   (rec["../data/" <> file, ##]& @@ opt),
   {file, {"coeff_prop_4l", "coeff_prop_4l_mod", "aajamp", "aajamp_mod"}},
   {opt, {{}, {"PolyRecMethod" -> FFPolyVandermonde}}}
];
