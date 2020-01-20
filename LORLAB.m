BeginPackage["LORLAB`"]

FFF[vars_, pars_][man_] := Module[{dumvars},
   dumvars = Table[Subscript[\[ScriptD], i], {i, 1, Length[vars]}];
   Transpose[D[man[dumvars, pars], {dumvars}]].D[
      man[dumvars, pars], {dumvars}] /. 
    Table[dumvars[[i]] -> vars[[i]], {i, 1, Length[vars]}]];

SFF[vars_, pars_][man_, normal_] := Module[{dumvars},
   dumvars = Table[Subscript[\[ScriptD], i], {i, 1, Length[vars]}];
   Table[normal[vars, pars].D[
       D[man[dumvars, pars], dumvars[[\[ScriptI]]]], 
       dumvars[[\[ScriptJ]]]], {\[ScriptI], 1, 
      Length[vars]}, {\[ScriptJ], 1, Length[vars]}] /. 
    Table[dumvars[[i]] -> vars[[i]], {i, 1, Length[vars]}]];

Zip[list1_, list2_] := If[Length[list1] != Length[list2],
   "List Mismatch",
   Sum[list1[[\[ScriptK]]]*list2[[\[ScriptK]]], {\[ScriptK], 1, 
     Length[list1]}]];

TE[vars_, pars_][man_, frame_] := 
  Module[{tanvars, normvars, dim, codim},
   codim = Length[frame];
   dim = Length[vars] - codim;
   tanvars = Take[vars, {1, dim}];
   normvars = Take[vars, {dim + 1, dim + codim}];
   FFF[tanvars, pars][man] + 
    Zip[normvars, SFF[tanvars, pars][man, #] & /@ frame]];

NE[vars_, pars_][man_, frame_] := 
  Module[{tanvars, dumTanVars, normvars, dim, codim},
   codim = Length[frame];
   dim = Length[vars] - codim;
   tanvars = Take[vars, {1, dim}];
   dumTanVars = Table[Subscript[\[ScriptD], i], {i, 1, dim}];
   normvars = Take[vars, {dim + 1, dim + codim}];
   Zip[normvars,
    	Table[#.D[normVec, {tanvars}], {normVec, #}] &[
      Table[vec[dumTanVars, pars], {vec, frame}]]
     	/. Table[
      dumTanVars[[\[ScriptK]]] -> tanvars[[\[ScriptK]]], {\[ScriptK], 
       1, dim}]]];
	
\[CapitalPsi][vars_, pars_][man_, frame_] := 
  Module[{tanVars, normVars, dim, coDim},
   coDim = Length[frame];
   dim = Length[vars] - coDim;
   tanVars = Take[vars, {1, dim}];
   normVars = Take[vars, {dim + 1, dim + coDim}];
   \[Sigma][tanVars, pars] + 
    Zip[normVars, Table[vec[tanVars, pars], {vec, frame}]]];

LORE[vars_, pars_][man_, frame_, field_] := 
  Module[{tanvars, dumTanVars, normvars, dim, codim},
   codim = Length[frame];
   dim = Length[vars] - codim;
   tanvars = Take[vars, {1, dim}];
   dumTanVars = Table[Subscript[\[ScriptD], i], {i, 1, dim}];
   normvars = Take[vars, {dim + 1, dim + codim}];
   Join[#1, 
       Table[vec[tanvars, pars], {vec, frame}].#2 - 
        NE[vars, pars][man, frame].#1] &[
     			Inverse[
       TE[vars, pars][man, 
        frame]].(field[\[CapitalPsi][vars, pars][man, frame], pars].D[
         man[dumTanVars, pars], {dumTanVars}]),
     			field[\[CapitalPsi][vars, pars][man, frame], pars]] /. 
    Table[dumTanVars[[\[ScriptK]]] -> 
      tanvars[[\[ScriptK]]], {\[ScriptK], 1, dim}]];
FLORE[vars_, pars_][vec_, basecurve_] := 
  Module[{w, nw, Nw, \[Kappa]w, list, n, K, Nvec},
   w[u_] := basecurve[u];
   nw[u_] := 
    Sqrt[D[basecurve[ry], ry].D[basecurve[ry], ry]] /. {ry -> u};
   n = Length[w[u]];
   list = FrenetSerretSystem[w[u], u];
   \[Kappa]w[p_, i_] := list[[1]][[i]] /. {u -> p};
   Nw[p_, i_] := list[[2]][[i]] /. {u -> p};
   K[p_] := 
    If[n == 2, 0, 
     DiagonalMatrix[Table[\[Kappa]w[p, i], {i, 2, n - 1}], 1] - 
      DiagonalMatrix[Table[\[Kappa]w[p, i], {i, 2, n - 1}], -1]];
   Nvec[vars2_, pars2_] := 
    Table[vec[
       w[vars2[[1]]] + Sum[vars2[[i]] Nw[vars2[[1]], i], {i, 2, n}], 
       pars2].Nw[vars2[[1]], j], {j, 2, n}];
   Join[{vec[
        w[vars[[1]]] + Sum[vars[[i]] Nw[vars[[1]], i], {i, 2, n}], 
        pars].Nw[vars[[1]], 
        1]/(nw[vars[[1]]] (1 - vars[[2]]*\[Kappa]w[vars[[1]], 1]))},
    Nvec[vars, 
      pars] + (vec[
          w[vars[[1]]] + Sum[vars[[i]] Nw[vars[[1]], i], {i, 2, n}], 
          pars].Nw[vars[[1]], 
          1]/((1 - vars[[2]]*\[Kappa]w[vars[[1]], 1]))) If[n == 2, 0, 
       K[vars[[1]]].Table[vars[[i]], {i, 2, n}]]]];

EndPackage[]