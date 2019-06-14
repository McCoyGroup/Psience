(* ::Package:: *)

pyReduce[expr_]:=
  Module[
    {e=expr},
    e=e/.{
        Rational[a_, b_]:>
          "("<>ToString[pyReduce[a]]<>")/("<>ToString[pyReduce[b]]<>")",
        Power[a_, b_?Negative]:>
          "1/("<>ToString[pyReduce[Power[a, Abs[b]]]]<>")",
        Power[a_, 1/2]:>
          "np.sqrt("<>ToString[pyReduce@a]<>")",
        HoldPattern[Power[a_, b_]]:>
          "("<>ToString[pyReduce@a]<>")"<>"**("<>ToString[pyReduce@b]<>")",
        HoldPattern[Times[-1, a__]]:>
          "-"<>pyReduce[Times[a]],
        HoldPattern[Times[a__]]:>StringRiffle[pyReduce/@{a},"*"],
        HoldPattern[Plus[a__]]:>StringRiffle[pyReduce/@{a},"+"],
        HoldPattern[Cos[a_]]:>
          "np.cos("<>ToString[pyReduce@a]<>")",
        HoldPattern[Sin[a_]]:>
          "np.sin("<>ToString[pyReduce@a]<>")",
        HoldPattern[ArcCos[a_]]:>
          "np.arccos("<>ToString[pyReduce@a]<>")",
        HoldPattern[ArcSin[a_]]:>
          "np.arcsin("<>ToString[pyReduce@a]<>")",
        HoldPattern[ArcTan[a_]]:>
          "np.arctan("<>ToString[pyReduce@a]<>")",
        HoldPattern[ArcTan[a_, b_]]:>
          "np.arctan2("<>ToString[pyReduce@a]<>", "<>ToString[pyReduce@b]<>")"
        };
    e
    ]


getPyEulerMat[orientation_]:=
  Module[
    {
      mat=EulerMatrix[{"0","1","2"}, orientation],
      strlens,
      paddedStrings
      },
    mat=mat/.{Cos[x_]:>"c["<>ToString[x]<>"]",Sin[x_]:>"s["<>ToString[x]<>"]"};
    mat=pyReduce[mat];
    strlens=Max/@Map[StringLength, Transpose[mat], {2}];
    paddedStrings=
      Transpose@
        MapThread[
          With[{l=#2}, Map[StringPadRight[#, l]&, #]]&,
          {
            Transpose[mat],
            strlens
            }
          ];
    "[\n"<>
    StringRiffle[
      Map["  [ "<>StringRiffle[#, ", "]<>" ]"&, paddedStrings],
      ",\n"
      ]<>"\n]"
    ]


getPyAnglesSystem//Clear
getPyAnglesSystem[orientation_]:=
  Module[
    {
      mat=EulerMatrix[{\[Alpha], \[Beta], \[Gamma]}, orientation],
      chobmat={{xx, xy, xz}, {yx, yy, yz}, {zx, zy, zz}},
      strlens,
      paddedStrings
      },
    Partition[Thread[Flatten[mat]==Flatten[chobmat]], 3]
    ]


formatPyAnglesSystem[solns_]:=
  Module[
    {
      baseSolns=Simplify@solns,
      redSolns
      },
    baseSolns=baseSolns/.
      Flatten[
        MapIndexed[
          #->"basis["<>ToString[#2[[1]]-1]<>"]["<>ToString[#2[[2]]-1]<>"]"&,
          {{xx, xy, xz}, {yx, yy, yz}, {zx, zy, zz}},
          {2}
          ]
        ];
    redSolns=pyReduce[baseSolns];
    "[\n"<>StringRiffle[redSolns, ",\n"]<>"\n]"
    ]


getPyAngles//Clear
getPyAngles[orientation_, 
  solEls_:{{3, 1}, {3, 3}, {1, 3}}, 
  format:True|False:True
  ]:=
  Module[
    {
      baseAngles=getPyAnglesSystem[orientation],
      baseSolns,
      noPiSolns
      },
    baseSolns=
      Solve[Extract[baseAngles, solEls], {\[Alpha], \[Beta], \[Gamma]}];
    baseSolns=
      Assuming[
        Flatten[{{xx, xy, xz}, {yx, yy, yz}, {zx, zy, zz}}]\[Element]Reals,
        {\[Alpha], \[Beta], \[Gamma]}/.baseSolns/._C->0//Simplify
        ]//Select[FreeQ[I|_Complex]];
    If[format,
      MinimalBy[
        formatPyAnglesSystem/@baseSolns,
        Max[StringLength[#]&/@#]&
        ][[1]], 
      baseSolns
      ]
    ]


(* ::Text:: *)
(*This one is somewhat subtle since you often need to choose your elements well. Therefore something like this is generally helpful as a preprocessing step:*)


getPyAnglesSystem[{1, 3, 2}]//TraditionalForm


getPyAngles[{1, 3, 2}, {{1, 2}, {2, 2}, {1, 1}}]


MinimalBy[getPyAngles[{1, 2, 3}, #]&/@Map[
Reverse,{
  {{1, 2}, {1, 3}, {3, 3}},
  {{1, 2}, {1, 1}, {3, 3}},
  {{1, 2}, {1, 1}, {2, 3}},
  {{1, 2}, {1, 3}, {2, 3}},
  {{1, 1}, {1, 3}, {2, 3}},
  {{1, 1}, {1, 3}, {3, 3}}
  },
  {2}
  ],
   StringLength
   ]


getPyEulerMat[{3, 2, 3}]
