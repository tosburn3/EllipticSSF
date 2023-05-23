(* ::Package:: *)

(* These first two definitions are prototypes for functions that will be overwritten 
when "InitializeSeff[]" is called from the linked C code for the effective source *)

PsiPm[m_, r_, \[Theta]_] := None; 

Seffm[m_, r_, \[Theta]_] := None; 


(* These global variables are prototypes that must be overwritten with actual 
values in the notebook where "EllipticSSF.wl" is imported via "Get[]" *)

r0 = None; 

a = None; 

v = None; (* v = 1/Sqrt[r0] *)

\[CapitalOmega] = None;  (* \[CapitalOmega] = v^3/(1+a v^3) *)


Delta[a_, rSubrPlus_] := Module[{r, rPlus, rMinus},
  rPlus = 1 + Sqrt[1 - a^2]; rMinus = 1 - Sqrt[1 - a^2]; r = rSubrPlus
     + rPlus; rSubrPlus (r - rMinus)
]


Sigma[a_, rSubrPlus_, \[Theta]_] := Module[{r, rPlus, rMinus},
  rPlus = 1 + Sqrt[1 - a^2]; rMinus = 1 - Sqrt[1 - a^2]; r = rSubrPlus
     + rPlus; Sqrt[(r^2 + a^2) ^ 2 - a^2 Delta[a, rSubrPlus] Sin[\[Theta]] ^ 2]
]


SigmaSq[a_, rSubrPlus_, \[Theta]_] := Module[{r, rPlus, rMinus},
  rPlus = 1 + Sqrt[1 - a^2]; rMinus = 1 - Sqrt[1 - a^2]; r = rSubrPlus
     + rPlus; (r^2 + a^2) ^ 2 - a^2 Delta[a, rSubrPlus] Sin[\[Theta]] ^ 2
]


\[CapitalDelta]\[Phi][a_, rSubrPlus_] := Module[{rPlus, rMinus, r},
  rPlus = 1 + Sqrt[1 - a^2]; rMinus = 1 - Sqrt[1 - a^2]; r = rSubrPlus
     + rPlus; a / (rPlus - rMinus) Log[rSubrPlus / (r - rMinus)]
]


d\[CapitalDelta]\[Phi]dr[a_, rSubrPlus_] := Module[{rPlus, rMinus, r},
  rPlus = 1 + Sqrt[1 - a^2]; rMinus = 1 - Sqrt[1 - a^2]; r = rSubrPlus
     + rPlus; a / (rSubrPlus (r - rMinus))
]


drStardr[a_, rSubrPlus_] := Module[{rPlus, rMinus, r},
  rPlus = 1 + Sqrt[1 - a^2]; r = rSubrPlus + rPlus; (r^2 + a^2) / Delta[
    a, rSubrPlus]
]


getrStarFromrSubrPlus[a_, rSubrPlus_] :=
  Module[{r, rPlus, rMinus},
    rPlus = 1 + Sqrt[1 - a^2]; rMinus = 1 - Sqrt[1 - a^2]; r = rSubrPlus
       + rPlus; r + 2 / (rPlus - rMinus) (rPlus Log[rSubrPlus / 2] - rMinus
       Log[(r - rMinus) / 2])
  ]


getrSubrPlusFromrStar[a_, rStar_] :=
  Module[{rSubrPlus, guess, rPlus},
    rPlus = 1 + Sqrt[1 - a^2];
    guess =
      If[rStar < -2,
        2 (1 - a^2) ^ (1 / rPlus - 1 / 2) Exp[(a^2 - rPlus + (rPlus -
           1) * rStar) / rPlus]
        ,
        If[rStar < 1000,
          2 * ProductLog[Exp[(rStar - rPlus) / 2]]
          ,
          rStar
        ]
      ];
    Re[
      rSubrPlus /.
        Quiet[
          Check[
            If[rStar < 10,
              FindRoot[getrStarFromrSubrPlus[a, rSubrPlus] - rStar, {
                rSubrPlus, guess}, AccuracyGoal -> 13, PrecisionGoal -> 13]
              ,
              FindRoot[getrStarFromrSubrPlus[a, rSubrPlus] / rStar - 
                1, {rSubrPlus, guess}, AccuracyGoal -> 13, PrecisionGoal -> 13]
            ]
            ,
            NSolve[getrStarFromrSubrPlus[a, rSubrPlus] / rStar - 1 ==
               0, rSubrPlus, Reals][[1]]
          ]
        ]
    ]
  ]


getrStarParams[a_, r0_, wtDiam_, rStarHguess_, rStarIguess_] :=
  Module[{rPlus, rStar0, rStarL, rStarR, rStarH, rStarI},
    rPlus = 1 + Sqrt[1 - a^2];
    rStar0 = getrStarFromrSubrPlus[a, r0 - rPlus];
    rStarL = rStar0 - 0.5 wtDiam;
    rStarR = rStar0 + 0.5 wtDiam;
    rStarI = rStarR + wtDiam Ceiling[(rStarIguess - rStarR) / wtDiam]
      ;
    rStarH = rStarL - wtDiam Ceiling[(rStarL - rStarHguess) / wtDiam]
      ;
    {rStarH, rStarL, rStar0, rStarR, rStarI}
  ]


getrStarListOld[rStarMax_, rStarMin_, \[CapitalDelta]rStar_] := Module[
  {},
 Table[rStar, {rStar, rStarMin, rStarMax, \[CapitalDelta]rStar}]
]


getrStarList[\[CapitalDelta]rStar_, rStarH_, rStarL_, rStar0_, rStarR_, rStarI_] :=
  Module[{wtDiam, wtNum, new\[CapitalDelta]rStar, totalNum},
    wtDiam = rStarR - rStarL;
    wtNum = 2 Round[(1 + wtDiam / \[CapitalDelta]rStar) / 2];
    new\[CapitalDelta]rStar = wtDiam / (wtNum - 1);
    totalNum = Round[(rStarI - rStarH) / new\[CapitalDelta]rStar] + 1;
    Return[Table[rStarH + (i - 1) new\[CapitalDelta]rStar, {i, 1, totalNum}]]
  ]


getrSubrPlusList[a_, rStarList_] := Table[getrSubrPlusFromrStar[a, rStar
  ], {rStar, rStarList}]


getIsourceBounds[rStarList_, rStarH_, rStarL_, rStar0_, rStarR_, rStarI_
  ] := Module[{\[CapitalDelta]rStar},
  \[CapitalDelta]rStar = (rStarList[[-1]] - rStarList[[1]]) / (Length[rStarList] - 
    1); Return[{1 + Round[(rStarL - rStarH) / \[CapitalDelta]rStar], 1 + Round[(rStarR 
    - rStarH) / \[CapitalDelta]rStar]}]
]


getThetaListOld[\[CapitalDelta]\[Theta]_] := Module[{},
Table[\[Theta], {\[Theta], 0, \[Pi], \[CapitalDelta]\[Theta]}]]


getThetaList[\[CapitalDelta]\[Theta]_, thetaSourceSize_] := Module[{numSourceInPi, newThetaSourceSize,
   sourceNum, nTheta, new\[CapitalDelta]\[Theta], totalNum},
  numSourceInPi = If[OddQ[Round[\[Pi] / thetaSourceSize]] == True,
    Round[\[Pi] / thetaSourceSize]
    ,
    2 (Round[\[Pi] / (2 thetaSourceSize)] + 1 / 2)
  ];
  newThetaSourceSize = N[\[Pi] / numSourceInPi];
  sourceNum = 2 Round[(1 + newThetaSourceSize / \[CapitalDelta]\[Theta]) / 2];
  new\[CapitalDelta]\[Theta] = newThetaSourceSize / (sourceNum - 1);
  totalNum = Round[\[Pi] / new\[CapitalDelta]\[Theta]] + 1;
  Return[Table[(j - 1) new\[CapitalDelta]\[Theta], {j, 1, totalNum}]]
]


getJsourceBounds[\[CapitalDelta]\[Theta]_, thetaSourceSize_] := Module[{numSourceInPi, newThetaSourceSize,
   sourceNum, nTheta, new\[CapitalDelta]\[Theta], totalNum},
  numSourceInPi = If[OddQ[Round[\[Pi] / thetaSourceSize]] == True,
    Round[\[Pi] / thetaSourceSize]
    ,
    2 (Round[\[Pi] / (2 thetaSourceSize)] + 1 / 2)
  ];
  newThetaSourceSize = N[\[Pi] / numSourceInPi];
  sourceNum = 2 Round[(1 + newThetaSourceSize / \[CapitalDelta]\[Theta]) / 2];
  new\[CapitalDelta]\[Theta] = newThetaSourceSize / (sourceNum - 1);
  totalNum = Round[\[Pi] / new\[CapitalDelta]\[Theta]] + 1;
  Return[{1 + Round[(\[Pi] / 2 - newThetaSourceSize / 2) / new\[CapitalDelta]\[Theta]], 1 + Round[
    (\[Pi] / 2 + newThetaSourceSize / 2) / new\[CapitalDelta]\[Theta]]}]
]


rowToI[row_, jMax_] := Floor[(row - 1) / jMax] + 1


rowToJ[row_, i_, jMax_] := row - (i - 1) jMax


ijToCol[i_, j_, jMax_] := j + (i - 1) jMax


fij[m_, a_, rSubrPlus_, \[Theta]_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, \[Omega]_] :=
  Module[{r, rPlus, \[CapitalSigma], \[CapitalDelta]},
    rPlus = 1 + Sqrt[1 - a^2];
    r = rSubrPlus + rPlus;
    \[CapitalSigma] = Sigma[a, rSubrPlus, \[Theta]];
    \[CapitalDelta] = Delta[a, rSubrPlus];
    -\[Omega]^2 + 4 a m \[Omega] r / \[CapitalSigma]^2 - (r^2 + a^2) ^ 2 / \[CapitalSigma]^2 (-2 / \[CapitalDelta]rStar^2) - 
      \[CapitalDelta] / \[CapitalSigma]^2 ((-2 / \[CapitalDelta]\[Theta]^2) - m^2 / Sin[\[Theta]] ^ 2 - 2 / r (1 - a^2 / r) - 2 I a
       m / r)
  ]


fiMj[m_, a_, rSubrPlus_, \[Theta]_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, \[Omega]_] :=
  Module[{r, rPlus, \[CapitalSigma], \[CapitalDelta]},
    rPlus = 1 + Sqrt[1 - a^2];
    r = rSubrPlus + rPlus;
    \[CapitalSigma] = Sigma[a, rSubrPlus, \[Theta]];
    \[CapitalDelta] = Delta[a, rSubrPlus];
    -(r^2 + a^2) ^ 2 / \[CapitalSigma]^2 (1 / \[CapitalDelta]rStar^2) - (2 I a m r (r^2 + a^2) - 
      2 a^2 \[CapitalDelta]) / (r \[CapitalSigma]^2) (-1 / (2 \[CapitalDelta]rStar))
  ]


fiPj[m_, a_, rSubrPlus_, \[Theta]_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, \[Omega]_] :=
  Module[{r, rPlus, \[CapitalSigma], \[CapitalDelta]},
    rPlus = 1 + Sqrt[1 - a^2];
    r = rSubrPlus + rPlus;
    \[CapitalSigma] = Sigma[a, rSubrPlus, \[Theta]];
    \[CapitalDelta] = Delta[a, rSubrPlus];
    -(r^2 + a^2) ^ 2 / \[CapitalSigma]^2 (1 / \[CapitalDelta]rStar^2) - (2 I a m r (r^2 + a^2) - 
      2 a^2 \[CapitalDelta]) / (r \[CapitalSigma]^2) (1 / (2 \[CapitalDelta]rStar))
  ]


fijM[m_, a_, rSubrPlus_, \[Theta]_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, \[Omega]_] :=
  Module[{r, rPlus, \[CapitalSigma], \[CapitalDelta]},
    rPlus = 1 + Sqrt[1 - a^2];
    r = rSubrPlus + rPlus;
    \[CapitalSigma] = Sigma[a, rSubrPlus, \[Theta]];
    \[CapitalDelta] = Delta[a, rSubrPlus];
    -\[CapitalDelta] / \[CapitalSigma]^2 ((1 / \[CapitalDelta]\[Theta]^2) + Cot[\[Theta]] (-1 / (2 \[CapitalDelta]\[Theta])))
  ]


fijP[m_, a_, rSubrPlus_, \[Theta]_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, \[Omega]_] :=
  Module[{r, rPlus, \[CapitalSigma], \[CapitalDelta]},
    rPlus = 1 + Sqrt[1 - a^2];
    r = rSubrPlus + rPlus;
    \[CapitalSigma] = Sigma[a, rSubrPlus, \[Theta]];
    \[CapitalDelta] = Delta[a, rSubrPlus];
    -\[CapitalDelta] / \[CapitalSigma]^2 ((1 / \[CapitalDelta]\[Theta]^2) + Cot[\[Theta]] (1 / (2 \[CapitalDelta]\[Theta])))
  ]


d2thetaCoeff[m_, a_, rSubrPlus_, \[Theta]_, \[Omega]_] :=
  -Delta[a, rSubrPlus] / SigmaSq[a, rSubrPlus, \[Theta]]


d1thetaCoeff[m_, a_, rSubrPlus_, \[Theta]_, \[Omega]_] :=
  - Delta[a, rSubrPlus] / SigmaSq[a, rSubrPlus, \[Theta]](Cot[\[Theta]])


d2rStarCoeff[m_, a_, rSubrPlus_, \[Theta]_, \[Omega]_] := Module[{r, rPlus},
    rPlus = 1 + Sqrt[1 - a^2];
    r = rSubrPlus + rPlus;
    -(r^2 + a^2) ^ 2 / SigmaSq[a, rSubrPlus, \[Theta]]]


d1rStarCoeff[m_, a_, rSubrPlus_, \[Theta]_, \[Omega]_] := Module[{r, rPlus},
    rPlus = 1 + Sqrt[1 - a^2];
    r = rSubrPlus + rPlus;
    - (2 I a m r (r^2 + a^2) - 
      2 a^2 Delta[a, rSubrPlus]) / (r SigmaSq[a, rSubrPlus, \[Theta]]) ]


d0Coeff[m_, a_, rSubrPlus_, \[Theta]_, \[Omega]_] :=
  Module[{r, rPlus},
    rPlus = 1 + Sqrt[1 - a^2];
    r = rSubrPlus + rPlus;
    -\[Omega]^2 + 4 a m \[Omega] r / SigmaSq[a, rSubrPlus, \[Theta]] - Delta[a, rSubrPlus]
       / SigmaSq[a, rSubrPlus, \[Theta]] (-m^2 / Sin[\[Theta]] ^ 2 - 2 / r (1 - a^2 / r) 
      - 2 I a m / r)
  ]


d1SecP1[h_] :=
  1 / (2 h)


d1SecM1[h_] :=
  -1 / (2 h)


d2SecP1[h_] :=
  1 / h^2


d2SecM1[h_] :=
  1 / h^2


d2Sec0[h_] :=
  -2 / h^2


d1FourthP1[h_] :=
  8 / (12 h)


d1FourthM1[h_] :=
  -8 / (12 h)


d1FourthP2[h_] :=
  -1 / (12 h)


d1FourthM2[h_] :=
  1 / (12 h)


d2FourthP1[h_] :=
  16 / (12 h^2)


d2FourthM1[h_] :=
  16 / (12 h^2)


d2FourthP2[h_] :=
  -1 / (12 h^2)


d2FourthM2[h_] :=
  -1 / (12 h^2)


d2Fourth0[h_] :=
  -30 / (12 h^2)


d1SkewM1[h_]:=-3/(12 h)


d1Skew0[h_]:=-10/(12 h)


d1SkewP1[h_]:=18/(12 h)


d1SkewP2[h_]:=-6/(12 h)


d1SkewP3[h_]:=1/(12 h)


d2SkewM1[h_]:=11/(12 h^2)


d2Skew0[h_]:=-20/(12 h^2)


d2SkewP1[h_]:=6/(12 h^2)


d2SkewP2[h_]:=4/(12 h^2)


d2SkewP3[h_]:=-1/(12 h^2)


couplingMatrixSlow[m_, rStarList_, rSubrPlusList_, thetaList_] :=
  Module[{iMax, jMax, i, j, mat, row, col, rSubrPlus, \[Theta], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], 
    \[Omega]},
    \[Omega] = m \[CapitalOmega];
    iMax = Length[rStarList];
    jMax = Length[thetaList];
    \[CapitalDelta]rStar = (rStarList[[-1]] - rStarList[[1]]) / (iMax - 1);
    \[CapitalDelta]\[Theta] = (thetaList[[-1]] - thetaList[[1]]) / (jMax - 1);
    mat = SparseArray[{}, {iMax jMax, iMax jMax}, 0];
    Do[
      \[Theta] = thetaList[[j]];
      Do[
        rSubrPlus = rSubrPlusList[[i]];
        row = ijToCol[i, j, jMax];
        mat[[row, row]] = fij[m, a, rSubrPlus, \[Theta], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], \[Omega]];
        col = ijToCol[i + 1, j, jMax];
        mat[[row, col]] = fiPj[m, a, rSubrPlus, \[Theta], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], \[Omega]];
        col = ijToCol[i - 1, j, jMax];
        mat[[row, col]] = fiMj[m, a, rSubrPlus, \[Theta], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], \[Omega]];
        col = ijToCol[i, j + 1, jMax];
        mat[[row, col]] = fijP[m, a, rSubrPlus, \[Theta], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], \[Omega]];
        col = ijToCol[i, j - 1, jMax];
        mat[[row, col]] = fijM[m, a, rSubrPlus, \[Theta], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], \[Omega]];
        
        ,
        {i, 2, iMax - 1}
      ];
      
      ,
      {j, 2, jMax - 1}
    ];
    If[m == 0,
      row = ijToCol[1, 1, jMax];
      mat[[row, row]] = 1;
      row = ijToCol[iMax, 1, jMax];
      mat[[row, row]] = 1;
      row = ijToCol[1, jMax, jMax];
      mat[[row, row]] = 1;
      row = ijToCol[iMax, jMax, jMax];
      mat[[row, row]] = 1;
      Do[
        row = ijToCol[i, 1, jMax];
        mat[[row, row]] = -3 / (2 \[CapitalDelta]\[Theta]);
        col = ijToCol[i, 2, jMax];
        mat[[row, col]] = 2 / \[CapitalDelta]\[Theta];
        col = ijToCol[i, 3, jMax];
        mat[[row, col]] = -1 / (2 \[CapitalDelta]\[Theta]);
        row = ijToCol[i, jMax, jMax];
        mat[[row, row]] = -3 / (2 \[CapitalDelta]\[Theta]);
        col = ijToCol[i, jMax - 1, jMax];
        mat[[row, col]] = 2 / \[CapitalDelta]\[Theta];
        col = ijToCol[i, jMax - 2, jMax];
        mat[[row, col]] = -1 / (2 \[CapitalDelta]\[Theta]);
        
        ,
        {i, 2, iMax - 1}
      ];
      Do[
        row = ijToCol[1, j, jMax];
        mat[[row, row]] = -3 / (2 \[CapitalDelta]rStar);
        col = ijToCol[2, j, jMax];
        mat[[row, col]] = 2 / \[CapitalDelta]rStar;
        col = ijToCol[3, j, jMax];
        mat[[row, col]] = -1 / (2 \[CapitalDelta]rStar);
        row = ijToCol[iMax, j, jMax];
        mat[[row, row]] = 1;
        
        ,
        {j, 2, jMax - 1}
      ];
      
      ,
      Do[
        row = ijToCol[i, 1, jMax];
        mat[[row, row]] = 1;
        row = ijToCol[i, jMax, jMax];
        mat[[row, row]] = 1;
        
        ,
        {i, 1, iMax}
      ];
      Do[
        row = ijToCol[1, j, jMax];
        mat[[row, row]] = I \[Omega] - 3 / (2 \[CapitalDelta]rStar);
        col = ijToCol[2, j, jMax];
        mat[[row, col]] = 2 / \[CapitalDelta]rStar;
        col = ijToCol[3, j, jMax];
        mat[[row, col]] = -1 / (2 \[CapitalDelta]rStar);
        row = ijToCol[iMax, j, jMax];
        mat[[row, row]] = -(I + \[CapitalDelta]rStar \[Omega]) (2 I + \[CapitalDelta]rStar \[Omega]) / \[CapitalDelta]rStar^2
          ;
        col = ijToCol[iMax - 1, j, jMax];
        mat[[row, col]] = (-5 + 4 I \[CapitalDelta]rStar \[Omega]) / \[CapitalDelta]rStar^2;
        col = ijToCol[iMax - 2, j, jMax];
        mat[[row, col]] = (4 - I \[CapitalDelta]rStar \[Omega]) / \[CapitalDelta]rStar^2;
        col = ijToCol[iMax - 3, j, jMax];
        mat[[row, col]] = -1 / \[CapitalDelta]rStar^2;
        
        ,
        {j, 2, jMax - 1}
      ];
      
    ];
    mat
  ]


couplingMatrix[m_, rStarList_, rSubrPlusList_, thetaList_] := Module[
  {iMax, jMax, i, j, mat, row, col, rSubrPlus, \[Theta], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], \[Omega]},
  \[Omega] = m \[CapitalOmega];
  iMax = Length[rStarList];
  jMax = Length[thetaList];
  \[CapitalDelta]rStar = (rStarList[[-1]] - rStarList[[1]]) / (iMax - 1);
  \[CapitalDelta]\[Theta] = (thetaList[[-1]] - thetaList[[1]]) / (jMax - 1);
  SparseArray[
    Join[
      If[m == 0,
        Join[{{1, 1} -> 1, {jMax, jMax} -> 1, Band[{2, 2}, {jMax - 1, jMax
           - 1}] -> -3 / (2 \[CapitalDelta]rStar), Band[{2, 2 + jMax}, {jMax - 1, 2 jMax - 1}]
           -> 2 / \[CapitalDelta]rStar, Band[{2, 2 + 2 jMax}, {jMax - 1, 3 jMax - 1}] -> -1 / (
          2 \[CapitalDelta]rStar), Band[{iMax jMax - 2 - (jMax - 3), iMax jMax - 2 - (jMax - 
          3)}, {iMax jMax, iMax jMax}] -> 1}, ParallelTable[{ijToCol[i, 1, jMax],
           ijToCol[i, 1, jMax]} -> -3 / (2 \[CapitalDelta]\[Theta]), {i, 2, iMax - 1}], ParallelTable[
          {ijToCol[i, 1, jMax], ijToCol[i, 2, jMax]} -> 2 / \[CapitalDelta]\[Theta], {i, 2, iMax - 1}
          ], ParallelTable[{ijToCol[i, 1, jMax], ijToCol[i, 3, jMax]} -> -1 / (2
           \[CapitalDelta]\[Theta]), {i, 2, iMax - 1}], ParallelTable[{ijToCol[i, jMax, jMax], ijToCol[
          i, jMax, jMax]} -> -3 / (2 \[CapitalDelta]\[Theta]), {i, 2, iMax - 1}], ParallelTable[{ijToCol[
          i, jMax, jMax], ijToCol[i, jMax - 1, jMax]} -> 2 / \[CapitalDelta]\[Theta], {i, 2, iMax - 1
          }], ParallelTable[{ijToCol[i, jMax, jMax], ijToCol[i, jMax - 2, jMax]
          } -> -1 / (2 \[CapitalDelta]\[Theta]), {i, 2, iMax - 1}]]
        ,
        Join[{Band[{2, 2}, {jMax - 1, jMax - 1}] -> I \[Omega] - 3 / (2 \[CapitalDelta]rStar
          ), Band[{2, ijToCol[2, 2, jMax]}, {jMax - 1, ijToCol[2, 2, jMax] + jMax
           - 3}] -> 2 / \[CapitalDelta]rStar, Band[{2, ijToCol[3, 2, jMax]}, {jMax - 1, ijToCol[
          3, 2, jMax] + jMax - 3}] -> -1 / (2 \[CapitalDelta]rStar), Band[{iMax jMax - 1 - (jMax
           - 3), iMax jMax - 1 - (jMax - 3)}, {iMax jMax - 1, iMax jMax - 1}] ->
           -(I + \[CapitalDelta]rStar \[Omega]) (2 I + \[CapitalDelta]rStar \[Omega]) / \[CapitalDelta]rStar^2, Band[{iMax jMax - 1 - (
          jMax - 3), ijToCol[iMax - 1, 2, jMax]}, {iMax jMax - 1, ijToCol[iMax 
          - 1, 2, jMax] + jMax - 3}] -> (-5 + 4 I \[CapitalDelta]rStar \[Omega]) / \[CapitalDelta]rStar^2, Band[{iMax
           jMax - 1 - (jMax - 3), ijToCol[iMax - 2, 2, jMax]}, {iMax jMax - 1, 
          ijToCol[iMax - 2, 2, jMax] + jMax - 3}] -> (4 - I \[CapitalDelta]rStar \[Omega]) / \[CapitalDelta]rStar^2,
           Band[{iMax jMax - 1 - (jMax - 3), ijToCol[iMax - 3, 2, jMax]}, {iMax
           jMax - 1, ijToCol[iMax - 3, 2, jMax] + jMax - 3}] -> -1 / \[CapitalDelta]rStar^2}, 
          ParallelTable[{ijToCol[i, 1, jMax], ijToCol[i, 1, jMax]} -> 1, {i, 1, 
          iMax}], ParallelTable[{ijToCol[i, jMax, jMax], ijToCol[i, jMax, jMax]
          } -> 1, {i, 1, iMax}]]
      ]
      ,
      ParallelTable[Band[{ijToCol[i, 2, jMax], ijToCol[i, 2, jMax]}] 
        -> fij[m, a, rSubrPlusList[[i]], thetaList[[2 ;; -2]], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], \[Omega]],
         {i, 2, iMax - 1}]
      ,
      ParallelTable[Band[{ijToCol[i, 2, jMax], ijToCol[i + 1, 2, jMax
        ]}] -> fiPj[m, a, rSubrPlusList[[i]], thetaList[[2 ;; -2]], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta],
         \[Omega]], {i, 2, iMax - 1}]
      ,
      ParallelTable[Band[{ijToCol[i, 2, jMax], ijToCol[i - 1, 2, jMax
        ]}] -> fiMj[m, a, rSubrPlusList[[i]], thetaList[[2 ;; -2]], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta],
         \[Omega]], {i, 2, iMax - 1}]
      ,
      ParallelTable[Band[{ijToCol[i, 2, jMax], ijToCol[i, 3, jMax]}] 
        -> fijP[m, a, rSubrPlusList[[i]], thetaList[[2 ;; -2]], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], \[Omega]],
         {i, 2, iMax - 1}]
      ,
      ParallelTable[Band[{ijToCol[i, 2, jMax], ijToCol[i, 1, jMax]}] 
        -> fijM[m, a, rSubrPlusList[[i]], thetaList[[2 ;; -2]], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], \[Omega]],
         {i, 2, iMax - 1}]
    ]
    ,
    {iMax jMax, iMax jMax}
    ,
    0
  ]
]


couplingMatrixFast[m_, rStarList_, rSubrPlusList_, thetaList_] :=
  Module[{M, iMax, jMax, i, j, mat, row, col, rSubrPlus, \[Theta], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta],
     \[Omega], rSubrPlusSlice, diag, rightDiag, leftDiag, rightSkipDiag, leftSkipDiag
    },
    M = 1;
    \[Omega] = m \[CapitalOmega];
    iMax = Length[rStarList];
    jMax = Length[thetaList];
    \[CapitalDelta]rStar = (rStarList[[-1]] - rStarList[[1]]) / (iMax - 1);
    \[CapitalDelta]\[Theta] = (thetaList[[-1]] - thetaList[[1]]) / (jMax - 1);
    rSubrPlusSlice = rSubrPlusList[[2 ;; -2]];
    diag = Flatten[ArrayPad[Transpose[Table[d0Coeff[m, a, rSubrPlusSlice,
       \[Theta], \[Omega]] + d2Sec0[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusSlice, \[Theta], \[Omega]] + d2Sec0[
      \[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusSlice, \[Theta], \[Omega]], {\[Theta], thetaList[[2 ;; -2]
      ]}]], 1]];
    rightDiag = Flatten[ArrayPad[Transpose[Table[d1SecP1[\[CapitalDelta]\[Theta]] d1thetaCoeff[
      m, a, rSubrPlusSlice, \[Theta], \[Omega]] + d2SecP1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusSlice,
       \[Theta], \[Omega]], {\[Theta], thetaList[[2 ;; -2]]}]], 1]];
    leftDiag = Flatten[ArrayPad[Transpose[Table[d1SecM1[\[CapitalDelta]\[Theta]] d1thetaCoeff[
      m, a, rSubrPlusSlice, \[Theta], \[Omega]] + d2SecM1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusSlice,
       \[Theta], \[Omega]], {\[Theta], thetaList[[2 ;; -2]]}]], 1]];
    rightSkipDiag = Flatten[ArrayPad[Transpose[Table[d1SecP1[\[CapitalDelta]rStar] 
      d1rStarCoeff[m, a, rSubrPlusSlice, \[Theta], \[Omega]] + d2SecP1[\[CapitalDelta]rStar] d2rStarCoeff[
      m, a, rSubrPlusSlice, \[Theta], \[Omega]], {\[Theta], thetaList[[2 ;; -2]]}]], {{1, 0},
       {1, 1}}]];
    leftSkipDiag = Flatten[ArrayPad[Transpose[Table[d1SecM1[\[CapitalDelta]rStar] 
      d1rStarCoeff[m, a, rSubrPlusSlice, \[Theta], \[Omega]] + d2SecM1[\[CapitalDelta]rStar] d2rStarCoeff[
      m, a, rSubrPlusSlice, \[Theta], \[Omega]], {\[Theta], thetaList[[2 ;; -2]]}]], {{
      0, 1}, {1, 1}}]];
    SparseArray[
      If[m == 0,
        Join[{{1, 1} -> 1, {jMax, jMax} -> 1, Band[{2, 2}, {jMax - 1,
           jMax - 1}] -> -3 / (2 \[CapitalDelta]rStar), Band[{2, 2 + jMax}, {jMax - 1, 2 jMax
           - 1}] -> 2 / \[CapitalDelta]rStar, Band[{2, 2 + 2 jMax}, {jMax - 1, 3 jMax - 1}] ->
           -1 / (2 \[CapitalDelta]rStar), Band[{iMax jMax - 2 - (jMax - 3), iMax jMax - 2 - (
          jMax - 3)}, {iMax jMax, iMax jMax}] -> 1}, Table[{ijToCol[i, 1, jMax],
           ijToCol[i, 1, jMax]} -> -3 / (2 \[CapitalDelta]\[Theta]), {i, 2, iMax - 1}], Table[{ijToCol[
          i, 1, jMax], ijToCol[i, 2, jMax]} -> 2 / \[CapitalDelta]\[Theta], {i, 2, iMax - 1}], Table[
          {ijToCol[i, 1, jMax], ijToCol[i, 3, jMax]} -> -1 / (2 \[CapitalDelta]\[Theta]), {i, 2, iMax
           - 1}], Table[{ijToCol[i, jMax, jMax], ijToCol[i, jMax, jMax]} -> -3 
          / (2 \[CapitalDelta]\[Theta]), {i, 2, iMax - 1}], Table[{ijToCol[i, jMax, jMax], ijToCol[i,
           jMax - 1, jMax]} -> 2 / \[CapitalDelta]\[Theta], {i, 2, iMax - 1}], Table[{ijToCol[i, jMax,
           jMax], ijToCol[i, jMax - 2, jMax]} -> -1 / (2 \[CapitalDelta]\[Theta]), {i, 2, iMax - 1}]
          ]
        ,
        Join[{Band[{2, 2}, {jMax - 1, jMax - 1}] -> I \[Omega] - 3 / (2 \[CapitalDelta]rStar
          ), Band[{2, ijToCol[2, 2, jMax]}, {jMax - 1, ijToCol[2, 2, jMax] + jMax
           - 3}] -> 2 / \[CapitalDelta]rStar, Band[{2, ijToCol[3, 2, jMax]}, {jMax - 1, ijToCol[
          3, 2, jMax] + jMax - 3}] -> -1 / (2 \[CapitalDelta]rStar), Band[{iMax jMax - 1 - (jMax
           - 3), iMax jMax - 1 - (jMax - 3)}, {iMax jMax - 1, iMax jMax - 1}] ->
           -(I + \[CapitalDelta]rStar \[Omega]) (2 I + \[CapitalDelta]rStar \[Omega]) / \[CapitalDelta]rStar^2, Band[{iMax jMax - 1 - (
          jMax - 3), ijToCol[iMax - 1, 2, jMax]}, {iMax jMax - 1, ijToCol[iMax 
          - 1, 2, jMax] + jMax - 3}] -> (-5 + 4 I \[CapitalDelta]rStar \[Omega]) / \[CapitalDelta]rStar^2, Band[{iMax
           jMax - 1 - (jMax - 3), ijToCol[iMax - 2, 2, jMax]}, {iMax jMax - 1, 
          ijToCol[iMax - 2, 2, jMax] + jMax - 3}] -> (4 - I \[CapitalDelta]rStar \[Omega]) / \[CapitalDelta]rStar^
          2, Band[{iMax jMax - 1 - (jMax - 3), ijToCol[iMax - 3, 2, jMax]}, {iMax
           jMax - 1, ijToCol[iMax - 3, 2, jMax] + jMax - 3}] -> -1 / \[CapitalDelta]rStar^2},
           Table[{ijToCol[i, 1, jMax], ijToCol[i, 1, jMax]} -> 1, {i, 1, iMax}],
           Table[{ijToCol[i, jMax, jMax], ijToCol[i, jMax, jMax]} -> 1, {i, 1, 
          iMax}]]
      ]
      ,
      {iMax jMax, iMax jMax}
      ,
      0
    ] + DiagonalMatrix[SparseArray[{i_} :> diag[[i]], {Length[diag]},
       0], 0] + DiagonalMatrix[SparseArray[{i_} :> rightDiag[[i]], {Length[
      rightDiag] - 1}, 0], 1] + DiagonalMatrix[SparseArray[{i_} :> leftDiag[[
      i + 1]], {Length[leftDiag] - 1}, 0], -1] + DiagonalMatrix[SparseArray[
      {i_} :> rightSkipDiag[[i]], {Length[rightSkipDiag]}, 0], jMax] + DiagonalMatrix[
      SparseArray[{i_} :> leftSkipDiag[[i]], {Length[leftSkipDiag]}, 0], -jMax
      ]
  ]


couplingMatrixFourth[m_, rStarList_, rSubrPlusList_, thetaList_] :=
  Module[{M, iMax, jMax, i, j, mat, row, col, rSubrPlus, \[Theta], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta],
     \[Omega], diag, rightDiag, leftDiag, rightSkipDiag, leftSkipDiag, rightDiag2,
     leftDiag2, rightSkipDiag2, leftSkipDiag2},
    M = 1;
    \[Omega] = m \[CapitalOmega];
    iMax = Length[rStarList];
    jMax = Length[thetaList];
    \[CapitalDelta]rStar = (rStarList[[-1]] - rStarList[[1]]) / (iMax - 1);
    \[CapitalDelta]\[Theta] = (thetaList[[-1]] - thetaList[[1]]) / (jMax - 1);
    diag = Flatten[ArrayPad[Transpose[Table[d0Coeff[m, a, rSubrPlusList[[
      3 ;; -3]], \[Theta], \[Omega]] + d2Fourth0[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusList[[
      3 ;; -3]], \[Theta], \[Omega]] + d2Fourth0[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[3 
      ;; -3]], \[Theta], \[Omega]], {\[Theta], thetaList[[3 ;; -3]]}]], 2]];
    rightDiag = Flatten[ArrayPad[Transpose[Table[d1FourthP1[\[CapitalDelta]\[Theta]] d1thetaCoeff[
      m, a, rSubrPlusList[[2 ;; -2]], \[Theta], \[Omega]] + d2FourthP1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m,
       a, rSubrPlusList[[2 ;; -2]], \[Theta], \[Omega]], {\[Theta], thetaList[[3 ;; -3]]}]], {{1,1},{2,2}}]]
      ;
    leftDiag = Flatten[ArrayPad[Transpose[Table[d1FourthM1[\[CapitalDelta]\[Theta]] d1thetaCoeff[
      m, a, rSubrPlusList[[2 ;; -2]], \[Theta], \[Omega]] + d2FourthM1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m,
       a, rSubrPlusList[[2 ;; -2]], \[Theta], \[Omega]], {\[Theta], thetaList[[3 ;; -3]]}]], {{1,1},{2,2}}]]
      ;
    rightDiag2 = Flatten[ArrayPad[Transpose[Table[d1FourthP2[\[CapitalDelta]\[Theta]] d1thetaCoeff[
      m, a, rSubrPlusList[[2 ;; -2]], \[Theta], \[Omega]] + d2FourthP2[\[CapitalDelta]\[Theta]] d2thetaCoeff[m,
       a, rSubrPlusList[[2 ;; -2]], \[Theta], \[Omega]], {\[Theta], thetaList[[3 ;; -3]]}]], {{1,1},{2,2}}]]
      ;
    leftDiag2 = Flatten[ArrayPad[Transpose[Table[d1FourthM2[\[CapitalDelta]\[Theta]] d1thetaCoeff[
      m, a, rSubrPlusList[[2 ;; -2]], \[Theta], \[Omega]] + d2FourthM2[\[CapitalDelta]\[Theta]] d2thetaCoeff[m,
       a, rSubrPlusList[[2 ;; -2]], \[Theta], \[Omega]], {\[Theta], thetaList[[3 ;; -3]]}]], {{1,1},{2,2}}]]
      ;
    rightSkipDiag = Flatten[ArrayPad[Transpose[Table[d1FourthP1[\[CapitalDelta]rStar
      ] d1rStarCoeff[m, a, rSubrPlusList[[3 ;; -3]], \[Theta], \[Omega]] + d2FourthP1[\[CapitalDelta]rStar
      ] d2rStarCoeff[m, a, rSubrPlusList[[3 ;; -3]], \[Theta], \[Omega]], {\[Theta], thetaList[[
      2 ;; -2]]}]], {{2, 1}, {1, 1}}]];
    leftSkipDiag = Flatten[ArrayPad[Transpose[Table[d1FourthM1[\[CapitalDelta]rStar
      ] d1rStarCoeff[m, a, rSubrPlusList[[3 ;; -3]], \[Theta], \[Omega]] + d2FourthM1[\[CapitalDelta]rStar
      ] d2rStarCoeff[m, a, rSubrPlusList[[3 ;; -3]], \[Theta], \[Omega]], {\[Theta], thetaList[[
      2 ;; -2]]}]], {{1, 2}, {1, 1}}]];
    rightSkipDiag2 = Flatten[ArrayPad[Transpose[Table[d1FourthP2[\[CapitalDelta]rStar
      ] d1rStarCoeff[m, a, rSubrPlusList[[3 ;; -3]], \[Theta], \[Omega]] + d2FourthP2[\[CapitalDelta]rStar
      ] d2rStarCoeff[m, a, rSubrPlusList[[3 ;; -3]], \[Theta], \[Omega]], {\[Theta], thetaList[[
      2 ;; -2]]}]], {{2, 0}, {1, 1}}]];
    leftSkipDiag2 = Flatten[ArrayPad[Transpose[Table[d1FourthM2[\[CapitalDelta]rStar
      ] d1rStarCoeff[m, a, rSubrPlusList[[3 ;; -3]], \[Theta], \[Omega]] + d2FourthM2[\[CapitalDelta]rStar
      ] d2rStarCoeff[m, a, rSubrPlusList[[3 ;; -3]], \[Theta], \[Omega]], {\[Theta], thetaList[[
      2 ;; -2]]}]], {{0, 2}, {1, 1}}]];
    SparseArray[
      Join[If[m == 0,
        Join[{{1, 1} -> 1, {jMax, jMax} -> 1, Band[{2, 2}, {jMax - 1,
           jMax - 1}] -> -25 / (12 \[CapitalDelta]rStar), Band[{2, 2 + jMax}, {jMax - 1, 2 jMax
           - 1}] -> 4 / \[CapitalDelta]rStar, Band[{2, 2 + 2 jMax}, {jMax - 1, 3 jMax - 1}] ->
           -3 / \[CapitalDelta]rStar, Band[{2, 2 + 3 jMax}, {jMax - 1, 4 jMax - 1}] -> 4 / (3
           \[CapitalDelta]rStar), Band[{2, 2 + 4 jMax}, {jMax - 1, 5 jMax - 1}] -> -1 / (4 \[CapitalDelta]rStar
          ), Band[{iMax jMax - 2 - (jMax - 3), iMax jMax - 2 - (jMax - 3)}, {iMax
           jMax, iMax jMax}] -> 1}, Table[{ijToCol[i, 1, jMax], ijToCol[i, 1, jMax
          ]} -> -25 / (12 \[CapitalDelta]\[Theta]), {i, 2, iMax - 1}], Table[{ijToCol[i, 1, jMax], ijToCol[
          i, 2, jMax]} -> 4 / \[CapitalDelta]\[Theta], {i, 2, iMax - 1}], Table[{ijToCol[i, 1, jMax],
           ijToCol[i, 3, jMax]} -> -3 / \[CapitalDelta]\[Theta], {i, 2, iMax - 1}], Table[{ijToCol[i,
           1, jMax], ijToCol[i, 4, jMax]} -> 4 / (3 \[CapitalDelta]\[Theta]), {i, 2, iMax - 1}], Table[
          {ijToCol[i, 1, jMax], ijToCol[i, 5, jMax]} -> -1 / (4 \[CapitalDelta]\[Theta]), {i, 2, iMax
           - 1}], Table[{ijToCol[i, jMax, jMax], ijToCol[i, jMax, jMax]} -> -25
           / (12 \[CapitalDelta]\[Theta]), {i, 2, iMax - 1}], Table[{ijToCol[i, jMax, jMax], ijToCol[
          i, jMax - 1, jMax]} -> 4 / \[CapitalDelta]\[Theta], {i, 2, iMax - 1}], Table[{ijToCol[i, jMax,
           jMax], ijToCol[i, jMax - 2, jMax]} -> -3 / \[CapitalDelta]\[Theta], {i, 2, iMax - 1}], Table[
          {ijToCol[i, jMax, jMax], ijToCol[i, jMax - 3, jMax]} -> 4 / (3 \[CapitalDelta]\[Theta]), {
          i, 2, iMax - 1}], Table[{ijToCol[i, jMax, jMax], ijToCol[i, jMax - 4,
           jMax]} -> -1 / (4 \[CapitalDelta]\[Theta]), {i, 2, iMax - 1}]]
        ,
        Join[{Band[{2, 2}, {jMax - 1, jMax - 1}] -> I \[Omega] - 25 / (12 \[CapitalDelta]rStar
          ), Band[{2, 2 + jMax}, {jMax - 1, 2 jMax - 1}] -> 4 / \[CapitalDelta]rStar, Band[{2,
           2 + 2 jMax}, {jMax - 1, 3 jMax - 1}] -> -3 / \[CapitalDelta]rStar, Band[{2, 2 + 3 
          jMax}, {jMax - 1, 4 jMax - 1}] -> 4 / (3 \[CapitalDelta]rStar), Band[{2, 2 + 4 jMax
          }, {jMax - 1, 5 jMax - 1}] -> -1 / (4 \[CapitalDelta]rStar), Band[{iMax jMax - 1 - 
          (jMax - 3), iMax jMax - 1 - (jMax - 3)}, {iMax jMax - 1, iMax jMax - 
          1}] -> \[CapitalDelta]rStar ^ (-4) - ((10 * I) * \[Omega]) / \[CapitalDelta]rStar^3 - (35 * \[Omega]^2) / (2 * 
          \[CapitalDelta]rStar^2) + (((25 * I) / 3) * \[Omega]^3) / \[CapitalDelta]rStar + \[Omega]^4, Band[{iMax jMax - 
          1 - (jMax - 3), ijToCol[iMax - 1, 2, jMax]}, {iMax jMax - 1, ijToCol[
          iMax - 1, 2, jMax] + jMax - 3}] -> -4 / \[CapitalDelta]rStar^4 + ((36 * I) * \[Omega]) / \[CapitalDelta]rStar
          ^3 + (52 * \[Omega]^2) / \[CapitalDelta]rStar^2 - ((16 * I) * \[Omega]^3) / \[CapitalDelta]rStar, Band[{iMax jMax
           - 1 - (jMax - 3), ijToCol[iMax - 2, 2, jMax]}, {iMax jMax - 1, ijToCol[
          iMax - 2, 2, jMax] + jMax - 3}] -> 6 / \[CapitalDelta]rStar^4 - ((48 * I) * \[Omega]) / \[CapitalDelta]rStar
          ^3 - (57 * \[Omega]^2) / \[CapitalDelta]rStar^2 + ((12 * I) * \[Omega]^3) / \[CapitalDelta]rStar, Band[{iMax jMax
           - 1 - (jMax - 3), ijToCol[iMax - 3, 2, jMax]}, {iMax jMax - 1, ijToCol[
          iMax - 3, 2, jMax] + jMax - 3}] -> -4 / \[CapitalDelta]rStar^4 + ((28 * I) * \[Omega]) / \[CapitalDelta]rStar
          ^3 + (28 * \[Omega]^2) / \[CapitalDelta]rStar^2 - (((16 * I) / 3) * \[Omega]^3) / \[CapitalDelta]rStar, Band[{iMax
           jMax - 1 - (jMax - 3), ijToCol[iMax - 4, 2, jMax]}, {iMax jMax - 1, 
          ijToCol[iMax - 4, 2, jMax] + jMax - 3}] -> \[CapitalDelta]rStar ^ (-4) - ((6 * I) *
           \[Omega]) / \[CapitalDelta]rStar^3 - (11 * \[Omega]^2) / (2 * \[CapitalDelta]rStar^2) + (I * \[Omega]^3) / \[CapitalDelta]rStar}, Table[
          {ijToCol[i, 1, jMax], ijToCol[i, 1, jMax]} -> 1, {i, 1, iMax}], Table[
          {ijToCol[i, jMax, jMax], ijToCol[i, jMax, jMax]} -> 1, {i, 1, iMax}]]
          
      ],
      {Band[{jMax + 2, 2}, {2 jMax-1, jMax-1}] -> d2SkewM1[\[CapitalDelta]rStar] 
      d2rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[2;;-2]], \[Omega]] + 
      d1SkewM1[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[2;;-2]], \[Omega]] ,
      Band[{jMax + 3, jMax + 3}, {2 jMax-2, 2 jMax-2}] -> d2Skew0[\[CapitalDelta]rStar] 
      d2rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[3;;-3]], \[Omega]] + 
      d1Skew0[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[3;;-3]], \[Omega]] +
      d0Coeff[m, a, rSubrPlusList[[2]], thetaList[[3;;-3]], \[Omega]] +
      d2Fourth0[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[2]], thetaList[[3;;-3]], \[Omega]] ,
      Band[{jMax + 2, 2 jMax + 2}, {2 jMax-1, 3 jMax-1}] -> d2SkewP1[\[CapitalDelta]rStar] 
      d2rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[2;;-2]], \[Omega]] + 
      d1SkewP1[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[2;;-2]], \[Omega]] ,
      Band[{jMax + 2, 3 jMax + 2}, {2 jMax-1, 4 jMax-1}] -> d2SkewP2[\[CapitalDelta]rStar] 
      d2rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[2;;-2]], \[Omega]] 
      + d1SkewP2[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[2;;-2]], \[Omega]] ,
      Band[{jMax + 2, 4 jMax + 2},{2 jMax-1, 5 jMax-1}] -> d2SkewP3[\[CapitalDelta]rStar] 
      d2rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[2;;-2]], \[Omega]] 
      + d1SkewP3[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[2;;-2]], \[Omega]] ,
      Band[{iMax jMax - 2 jMax+2,  iMax jMax - jMax+2},{iMax jMax - jMax-1, iMax jMax - 1}] ->
      d2SkewM1[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[2;;-2]], \[Omega]] - 
      d1SkewM1[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[2;;-2]], \[Omega]] ,
      Band[{iMax jMax - 2 jMax+3,  iMax jMax - 2 jMax+3},{iMax jMax - jMax-2, iMax jMax - jMax-2}] -> 
      d2Skew0[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[3;;-3]], \[Omega]] - 
      d1Skew0[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[3;;-3]], \[Omega]] +
      d0Coeff[m, a, rSubrPlusList[[-2]], thetaList[[3;;-3]], \[Omega]] +
      d2Fourth0[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[-2]], thetaList[[3;;-3]], \[Omega]] , 
      Band[{iMax jMax - 2 jMax+2,  iMax jMax - 3 jMax+2},{iMax jMax - jMax-1, iMax jMax - 2 jMax-1}] -> 
      d2SkewP1[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[2;;-2]], \[Omega]] - 
      d1SkewP1[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[2;;-2]], \[Omega]] ,
      Band[{iMax jMax - 2 jMax+2,  iMax jMax - 4 jMax+2},{iMax jMax - jMax-1, iMax jMax - 3 jMax-1}] -> 
      d2SkewP2[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[2;;-2]], \[Omega]] - 
      d1SkewP2[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[2;;-2]], \[Omega]] ,
      Band[{iMax jMax - 2 jMax+2,  iMax jMax - 5 jMax+2},{iMax jMax - jMax-1, iMax jMax - 4 jMax-1}] -> 
      d2SkewP3[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[2;;-2]], \[Omega]] - 
      d1SkewP3[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[2;;-2]], \[Omega]] ,
      {jMax + 2, jMax + 2} -> d0Coeff[m, a, rSubrPlusList[[2]], thetaList[[2]], \[Omega]] + 
      d2Skew0[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[2]], \[Omega]] + 
      d1Skew0[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[2]], \[Omega]] + 
      d2Skew0[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[2]], thetaList[[2]], \[Omega]] + 
      d1Skew0[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[2]], thetaList[[2]], \[Omega]] ,
      {2 jMax - 1, 2 jMax - 1} -> d0Coeff[m, a, rSubrPlusList[[2]], thetaList[[-2]], \[Omega]] + 
      d2Skew0[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[-2]], \[Omega]] + 
      d1Skew0[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[2]], thetaList[[-2]], \[Omega]] + 
      d2Skew0[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[2]], thetaList[[-2]], \[Omega]] - 
      d1Skew0[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[2]], thetaList[[-2]], \[Omega]] ,
      {iMax jMax - jMax-1, iMax jMax - jMax-1} -> d0Coeff[m, a, rSubrPlusList[[-2]], thetaList[[-2]], \[Omega]] + 
      d2Skew0[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[-2]], \[Omega]] - 
      d1Skew0[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[-2]], \[Omega]] + 
      d2Skew0[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[-2]], thetaList[[-2]], \[Omega]] - 
      d1Skew0[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[-2]], thetaList[[-2]], \[Omega]] ,
      {iMax jMax - 2 jMax+2, iMax jMax - 2 jMax+2} -> d0Coeff[m, a, rSubrPlusList[[-2]], thetaList[[2]], \[Omega]] + 
      d2Skew0[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[2]], \[Omega]] - 
      d1Skew0[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlusList[[-2]], thetaList[[2]], \[Omega]] + 
      d2Skew0[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[-2]], thetaList[[2]], \[Omega]] + 
      d1Skew0[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[-2]], thetaList[[2]], \[Omega]]
      },
      Table[{(i-1)jMax+2, (i-1)jMax+1} -> 
      d2SkewM1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[2]], \[Omega]] + 
      d1SkewM1[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[2]], \[Omega]], {i, 2, iMax-1}],
      Table[{(i-1)jMax+2, (i-1)jMax+2} -> d0Coeff[m, a, rSubrPlusList[[i]], thetaList[[2]], \[Omega]]+
      d2Fourth0[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusList[[i]], thetaList[[2]], \[Omega]]+
      d2Skew0[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[2]], \[Omega]] + 
      d1Skew0[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[2]], \[Omega]], {i, 3, iMax-2}],
      Table[{(i-1)jMax+2, (i-1)jMax+3} -> 
      d2SkewP1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[2]], \[Omega]] + 
      d1SkewP1[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[2]], \[Omega]], {i, 2, iMax-1}],
      Table[{(i-1)jMax+2, (i-1)jMax+4} -> 
      d2SkewP2[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[2]], \[Omega]] + 
      d1SkewP2[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[2]], \[Omega]], {i, 2, iMax-1}],
      Table[{(i-1)jMax+2, (i-1)jMax+5} -> 
      d2SkewP3[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[2]], \[Omega]] + 
      d1SkewP3[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[2]], \[Omega]], {i, 2, iMax-1}],
      Table[{i jMax-1, i jMax} -> 
      d2SkewM1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[-2]], \[Omega]] - 
      d1SkewM1[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[-2]], \[Omega]], {i, 2, iMax-1}],
      Table[{i jMax-1, i jMax-1} -> d0Coeff[m, a, rSubrPlusList[[i]], thetaList[[-2]], \[Omega]]+
      d2Fourth0[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlusList[[i]], thetaList[[-2]], \[Omega]]+
      d2Skew0[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[-2]], \[Omega]] - 
      d1Skew0[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[-2]], \[Omega]], {i, 3, iMax-2}],
      Table[{i jMax-1, i jMax-2} -> 
      d2SkewP1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[-2]], \[Omega]] - 
      d1SkewP1[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[-2]], \[Omega]], {i, 2, iMax-1}],
      Table[{i jMax-1, i jMax-3} -> 
      d2SkewP2[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[-2]], \[Omega]] - 
      d1SkewP2[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[-2]], \[Omega]], {i, 2, iMax-1}],
      Table[{i jMax-1, i jMax-4} -> 
      d2SkewP3[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[-2]], \[Omega]] - 
      d1SkewP3[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlusList[[i]], thetaList[[-2]], \[Omega]], {i, 2, iMax-1}]
      ]
      ,
      {iMax jMax, iMax jMax}
      ,
      0
    ]
     + DiagonalMatrix[SparseArray[{i_} :> diag[[i]], {Length[diag]},
       0], 0] + DiagonalMatrix[SparseArray[{i_} :> rightDiag[[i]], {Length[
      rightDiag] - 1}, 0], 1] + DiagonalMatrix[SparseArray[{i_} :> leftDiag[[
      i + 1]], {Length[leftDiag] - 1}, 0], -1] + DiagonalMatrix[SparseArray[
      {i_} :> rightDiag2[[i]], {Length[rightDiag2] - 2}, 0], 2] + DiagonalMatrix[
      SparseArray[{i_} :> leftDiag2[[i + 2]], {Length[leftDiag2] - 2}, 0], 
      -2] + DiagonalMatrix[SparseArray[{i_} :> rightSkipDiag[[i]], {Length[
      rightSkipDiag]}, 0], jMax] + DiagonalMatrix[SparseArray[{i_} :> leftSkipDiag[[
      i]], {Length[leftSkipDiag]}, 0], -jMax] + DiagonalMatrix[SparseArray[
      {i_} :> rightSkipDiag2[[i]], {Length[rightSkipDiag2]}, 0], 2 * jMax] 
      + DiagonalMatrix[SparseArray[{i_} :> leftSkipDiag2[[i]], {Length[leftSkipDiag2
      ]}, 0], -2 * jMax]
  ]


testMatrixSlow[m_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, rStarMax_, rStarMin_] := Module[{rStarList,
   rSubrPlusList, thetaList, mat},
  rStarList = getrStarListOld[rStarMax, rStarMin, \[CapitalDelta]rStar]; rSubrPlusList
     = getrSubrPlusList[a, rStarList]; thetaList = getThetaListOld[\[CapitalDelta]\[Theta]]; mat = couplingMatrixSlow[
    m, rStarList, rSubrPlusList, thetaList]; mat
]


testMatrix[m_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, rStarMax_, rStarMin_] := Module[{rStarList,
   rSubrPlusList, thetaList, mat},
  rStarList = getrStarListOld[rStarMax, rStarMin, \[CapitalDelta]rStar]; rSubrPlusList
     = getrSubrPlusList[a, rStarList]; thetaList = getThetaListOld[\[CapitalDelta]\[Theta]]; mat = couplingMatrix[
    m, rStarList, rSubrPlusList, thetaList]; mat
]


testMatrixFast[m_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, rStarMax_, rStarMin_] := Module[{rStarList,
   rSubrPlusList, thetaList, mat},
  rStarList = getrStarListOld[rStarMax, rStarMin, \[CapitalDelta]rStar]; rSubrPlusList
     = getrSubrPlusList[a, rStarList];
       thetaList = getThetaListOld[\[CapitalDelta]\[Theta]]; mat = couplingMatrixFast[
    m, rStarList, rSubrPlusList, thetaList]; mat
]


testMatrixFourth[m_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, rStarMax_, rStarMin_] := Module[{rStarList,
   rSubrPlusList, thetaList, mat},
  rStarList = getrStarListOld[rStarMax, rStarMin, \[CapitalDelta]rStar]; rSubrPlusList
     = getrSubrPlusList[a, rStarList];
       thetaList = getThetaListOld[\[CapitalDelta]\[Theta]]; mat = couplingMatrixFourth[
    m, rStarList, rSubrPlusList, thetaList]; mat
]


sourceVector[m_, rStarList_, rSubrPlusList_, thetaList_, iSourceMin_,
   iSourceMax_, jSourceMin_, jSourceMax_] :=
  Module[{iMax, jMax, i, j, vec, row, rSubrPlus, \[Theta], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], \[Omega], rPlus,
     r, fact, weight},
    weight = 1 / (2 \[Pi]);
    \[Omega] = m \[CapitalOmega];
    rPlus = 1 + Sqrt[1 - a^2];
    iMax = Length[rStarList];
    jMax = Length[thetaList];
    \[CapitalDelta]rStar = (rStarList[[-1]] - rStarList[[1]]) / (iMax - 1);
    \[CapitalDelta]\[Theta] = (thetaList[[-1]] - thetaList[[1]]) / (jMax - 1);
    vec = SparseArray[{}, {iMax jMax}, 0];
    Do[
      \[Theta] = thetaList[[j]];
      Do[
        vec[[ijToCol[i, j, jMax]]] = Exp[-I m \[CapitalDelta]\[Phi][a, rSubrPlusList[[i]
          ]]] Seffm[m, rSubrPlusList[[i]] + rPlus, \[Theta]] weight;
        
        ,
        {i, iSourceMin, iSourceMax}
      ];
      
      ,
      {j, jSourceMin, jSourceMax}
    ];
    Do[
      rSubrPlus = rSubrPlusList[[i]];
      fact = Exp[-I m \[CapitalDelta]\[Phi][a, rSubrPlus]];
      r = rSubrPlusList[[i]] + rPlus;
      \[Theta] = thetaList[[jSourceMin]];
      row = ijToCol[i, jSourceMin, jMax];
      vec[[row]] = vec[[row]] + (d1SecM1[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2SecM1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMin - 1]]] weight;
      \[Theta] = thetaList[[jSourceMax]];
      row = ijToCol[i, jSourceMax, jMax];
      vec[[row]] = vec[[row]] + (d1SecP1[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2SecP1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMax + 1]]] weight;
      \[Theta] = thetaList[[jSourceMin - 1]];
      row = ijToCol[i, jSourceMin - 1, jMax];
      vec[[row]] = vec[[row]] - (d1SecP1[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2SecP1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMin]]] weight;
      \[Theta] = thetaList[[jSourceMax + 1]];
      row = ijToCol[i, jSourceMax + 1, jMax];
      vec[[row]] = vec[[row]] - (d1SecM1[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2SecM1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMax]]] weight;
      
      ,
      {i, iSourceMin, iSourceMax}
    ];
    Do[
      \[Theta] = thetaList[[j]];
      rSubrPlus = rSubrPlusList[[iSourceMin]];
      row = ijToCol[iSourceMin, j, jMax];
      r = rSubrPlusList[[iSourceMin - 1]] + rPlus;
      fact = Exp[-I m \[CapitalDelta]\[Phi][a, r-rPlus]];
      vec[[row]] = vec[[row]] + (d1SecM1[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2SecM1[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, \[Theta]] weight;
      rSubrPlus = rSubrPlusList[[iSourceMax]];
      row = ijToCol[iSourceMax, j, jMax];
      r = rSubrPlusList[[iSourceMax + 1]] + rPlus;
      fact = Exp[-I m \[CapitalDelta]\[Phi][a, r-rPlus]];
      vec[[row]] = vec[[row]] + (d1SecP1[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2SecP1[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, \[Theta]] weight;
      rSubrPlus = rSubrPlusList[[iSourceMin - 1]];
      row = ijToCol[iSourceMin - 1, j, jMax];
      r = rSubrPlusList[[iSourceMin]] + rPlus;
      fact = Exp[-I m \[CapitalDelta]\[Phi][a, r-rPlus]];
      vec[[row]] = vec[[row]] - (d1SecP1[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2SecP1[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, \[Theta]] weight;
      rSubrPlus = rSubrPlusList[[iSourceMax + 1]];
      row = ijToCol[iSourceMax + 1, j, jMax];
      r = rSubrPlusList[[iSourceMax]] + rPlus;
      fact = Exp[-I m \[CapitalDelta]\[Phi][a, r-rPlus]];
      vec[[row]] = vec[[row]] - (d1SecM1[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2SecM1[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, \[Theta]] weight;
      
      ,
      {j, jSourceMin, jSourceMax}
    ];
    vec
  ]


sourceVectorFourth[m_, rStarList_, rSubrPlusList_, thetaList_, iSourceMin_,
   iSourceMax_, jSourceMin_, jSourceMax_] :=
  Module[{iMax, jMax, i, j, vec, row, rSubrPlus, \[Theta], \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], \[Omega], rPlus,
     r, fact, weight, r2, fact2},
    weight = 1 / (2 \[Pi]);
    \[Omega] = m \[CapitalOmega];
    rPlus = 1 + Sqrt[1 - a^2];
    iMax = Length[rStarList];
    jMax = Length[thetaList];
    \[CapitalDelta]rStar = (rStarList[[-1]] - rStarList[[1]]) / (iMax - 1);
    \[CapitalDelta]\[Theta] = (thetaList[[-1]] - thetaList[[1]]) / (jMax - 1);
    vec = SparseArray[{}, {iMax jMax}, 0];
    Do[
      \[Theta] = thetaList[[j]];
      Do[
        vec[[ijToCol[i, j, jMax]]] = Exp[-I m \[CapitalDelta]\[Phi][a, rSubrPlusList[[i]
          ]]] Seffm[m, rSubrPlusList[[i]] + rPlus, \[Theta]] weight;
        
        ,
        {i, iSourceMin, iSourceMax}
      ];
      
      ,
      {j, jSourceMin, jSourceMax}
    ];
    Do[
      rSubrPlus = rSubrPlusList[[i]];
      fact = Exp[-I m \[CapitalDelta]\[Phi][a, rSubrPlus]];
      r = rSubrPlusList[[i]] + rPlus;
      \[Theta] = thetaList[[jSourceMin]];
      row = ijToCol[i, jSourceMin, jMax];
      vec[[row]] = vec[[row]] + (d1FourthM1[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthM1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMin - 1]]] weight + 
        (d1FourthM2[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthM2[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMin - 2]]] weight;
      \[Theta] = thetaList[[jSourceMin+1]];
      row = ijToCol[i, jSourceMin+1, jMax];
      vec[[row]] = vec[[row]] + (d1FourthM2[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthM2[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMin - 1]]] weight;
      \[Theta] = thetaList[[jSourceMax]];
      row = ijToCol[i, jSourceMax, jMax];
      vec[[row]] = vec[[row]] + (d1FourthP1[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthP1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMax + 1]]] weight +
         (d1FourthP2[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthP2[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMax + 2]]] weight;
      \[Theta] = thetaList[[jSourceMax-1]];
      row = ijToCol[i, jSourceMax-1, jMax];
       vec[[row]] = vec[[row]] + (d1FourthP2[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthP2[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMax + 1]]] weight;
      \[Theta] = thetaList[[jSourceMin - 1]];
      row = ijToCol[i, jSourceMin - 1, jMax];
      vec[[row]] = vec[[row]] - (d1FourthP1[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthP1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMin]]] weight - 
        (d1FourthP2[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthP2[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMin+1]]] weight;
      \[Theta] = thetaList[[jSourceMin - 2]];
      row = ijToCol[i, jSourceMin - 2, jMax];
      vec[[row]] = vec[[row]] - (d1FourthP2[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthP2[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMin]]] weight;
      \[Theta] = thetaList[[jSourceMax + 1]];
      row = ijToCol[i, jSourceMax + 1, jMax];
      vec[[row]] = vec[[row]] - (d1FourthM1[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthM1[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMax]]] weight - 
        (d1FourthM2[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthM2[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMax-1]]] weight;
      \[Theta] = thetaList[[jSourceMax + 2]];
      row = ijToCol[i, jSourceMax + 2, jMax];
      vec[[row]] = vec[[row]] - (d1FourthM2[\[CapitalDelta]\[Theta]] d1thetaCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthM2[\[CapitalDelta]\[Theta]] d2thetaCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, thetaList[[jSourceMax]]] weight;
      
      ,
      {i, iSourceMin, iSourceMax}
    ];
    Do[
      \[Theta] = thetaList[[j]];
      rSubrPlus = rSubrPlusList[[iSourceMin]];
      row = ijToCol[iSourceMin, j, jMax];
      r = rSubrPlusList[[iSourceMin - 1]] + rPlus;
      r2 = rSubrPlusList[[iSourceMin - 2]] + rPlus;
      fact = Exp[-I m \[CapitalDelta]\[Phi][a, r-rPlus]];
      fact2 = Exp[-I m \[CapitalDelta]\[Phi][a, r2-rPlus]];
      vec[[row]] = vec[[row]] + (d1FourthM1[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthM1[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, \[Theta]] weight + (d1FourthM2[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthM2[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact2 PsiPm[
        m, r2, \[Theta]] weight;
      rSubrPlus = rSubrPlusList[[iSourceMin+1]];
      row = ijToCol[iSourceMin+1, j, jMax];
      r2 = rSubrPlusList[[iSourceMin - 1]] + rPlus;
      fact2 = Exp[-I m \[CapitalDelta]\[Phi][a, r2-rPlus]];
      vec[[row]] = vec[[row]] + (d1FourthM2[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthM2[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact2 PsiPm[
        m, r2, \[Theta]] weight;
      rSubrPlus = rSubrPlusList[[iSourceMax]];
      row = ijToCol[iSourceMax, j, jMax];
      r = rSubrPlusList[[iSourceMax + 1]] + rPlus;
      r2 = rSubrPlusList[[iSourceMax + 2]] + rPlus;
      fact = Exp[-I m \[CapitalDelta]\[Phi][a, r-rPlus]];
      fact2 = Exp[-I m \[CapitalDelta]\[Phi][a, r2-rPlus]];
      vec[[row]] = vec[[row]] + (d1FourthP1[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthP1[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, \[Theta]] weight + (d1FourthP2[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthP2[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact2 PsiPm[
        m, r2, \[Theta]] weight;
      rSubrPlus = rSubrPlusList[[iSourceMax-1]];
      row = ijToCol[iSourceMax-1, j, jMax];
      r2 = rSubrPlusList[[iSourceMax + 1]] + rPlus;
      fact2 = Exp[-I m \[CapitalDelta]\[Phi][a, r2-rPlus]];
      vec[[row]] = vec[[row]] + (d1FourthP2[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthP2[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact2 PsiPm[
        m, r2, \[Theta]] weight;
      rSubrPlus = rSubrPlusList[[iSourceMin - 1]];
      row = ijToCol[iSourceMin - 1, j, jMax];
      r = rSubrPlusList[[iSourceMin]] + rPlus;
      r2 = rSubrPlusList[[iSourceMin + 1]] + rPlus;
      fact = Exp[-I m \[CapitalDelta]\[Phi][a, r-rPlus]];
      fact2 = Exp[-I m \[CapitalDelta]\[Phi][a, r2-rPlus]];
      vec[[row]] = vec[[row]] - (d1FourthP1[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthP1[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, \[Theta]] weight - (d1FourthP2[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthP2[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact2 PsiPm[
        m, r2, \[Theta]] weight;
      rSubrPlus = rSubrPlusList[[iSourceMin - 2]];
      row = ijToCol[iSourceMin - 2, j, jMax];
      r2 = rSubrPlusList[[iSourceMin]] + rPlus;
      fact2 = Exp[-I m \[CapitalDelta]\[Phi][a, r2-rPlus]];
      vec[[row]] = vec[[row]] - (d1FourthP2[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthP2[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact2 PsiPm[
        m, r2, \[Theta]] weight;
      rSubrPlus = rSubrPlusList[[iSourceMax + 1]];
      row = ijToCol[iSourceMax + 1, j, jMax];
      r = rSubrPlusList[[iSourceMax]] + rPlus;
      r2 = rSubrPlusList[[iSourceMax - 1]] + rPlus;
      fact = Exp[-I m \[CapitalDelta]\[Phi][a, r-rPlus]];
      fact2 = Exp[-I m \[CapitalDelta]\[Phi][a, r2-rPlus]];
      vec[[row]] = vec[[row]] - (d1FourthM1[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthM1[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact PsiPm[
        m, r, \[Theta]] weight - (d1FourthM2[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthM2[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact2 PsiPm[
        m, r2, \[Theta]] weight;
      rSubrPlus = rSubrPlusList[[iSourceMax + 2]];
      row = ijToCol[iSourceMax + 2, j, jMax];
      r2 = rSubrPlusList[[iSourceMax]] + rPlus;
      fact2 = Exp[-I m \[CapitalDelta]\[Phi][a, r2-rPlus]];
      vec[[row]] = vec[[row]] - (d1FourthM2[\[CapitalDelta]rStar] d1rStarCoeff[m, a, rSubrPlus,
         \[Theta], \[Omega]] + d2FourthM2[\[CapitalDelta]rStar] d2rStarCoeff[m, a, rSubrPlus, \[Theta], \[Omega]]) fact2 PsiPm[
        m, r2, \[Theta]] weight;
      
      ,
      {j, jSourceMin, jSourceMax}
    ];
    vec
  ]


testSourceVector[m_, wtDiam_, thetaSourceSize_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, rStarHguess_,
   rStarIguess_] :=
  Module[{rStarList, rSubrPlusList, thetaList, iSourceMin, iSourceMax,
     jSourceMin, jSourceMax, vec, rStarH, rStarL, rStar0, rStarR, rStarI,
     jMax, source2D, Phi2D, rPlus},
    {rStarH, rStarL, rStar0, rStarR, rStarI} = getrStarParams[a, r0, 
      wtDiam, rStarHguess, rStarIguess];
    rStarList = getrStarList[\[CapitalDelta]rStar, rStarH, rStarL, rStar0, rStarR, 
      rStarI];
    rSubrPlusList = getrSubrPlusList[a, rStarList];
    thetaList = getThetaList[\[CapitalDelta]\[Theta], thetaSourceSize];
    {iSourceMin, iSourceMax} = getIsourceBounds[rStarList, rStarH, rStarL,
       rStar0, rStarR, rStarI];
    {jSourceMin, jSourceMax} = getJsourceBounds[\[CapitalDelta]\[Theta], thetaSourceSize];
      
    vec = sourceVector[m, rStarList, rSubrPlusList, thetaList, iSourceMin,
       iSourceMax, jSourceMin, jSourceMax];
    jMax = Length[thetaList];
    source2D = Partition[vec, jMax][[iSourceMin + 1 ;; iSourceMax - 1,
       jSourceMin + 1 ;; jSourceMax - 1]];
    rPlus = 1 + Sqrt[1 - a^2];
    Phi2D = Table[Exp[-I m \[CapitalDelta]\[Phi][a, rSubrPlusList[[i]]]] PsiPm[m, rSubrPlusList[[
      i]] + rPlus, thetaList[[j]]] / (2 \[Pi]), {i, iSourceMin + 1, iSourceMax 
      - 1}, {j, jSourceMin + 1, jSourceMax - 1}];
    Return[{source2D, Phi2D, rStarList, rSubrPlusList, thetaList, iSourceMin,
       iSourceMax, jSourceMin, jSourceMax}]
  ]


getFm[sol_, m_, rStarList_, rSubrPlusList_, thetaList_, iSourceMin_, 
  iSourceMax_, jSourceMin_, jSourceMax_] :=
  Module[{solSource, \[CapitalPsi], d\[CapitalPsi]drStar, d\[CapitalPsi]d\[Theta], interp, interpDrStar, interpD\[Theta],
     rStar, \[Theta], \[Omega], rPlus, iMax, jMax, \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], rStar0, r0, rSubrPlus0, 
    d\[CapitalPhi]dr, d\[CapitalPhi]dt, d\[CapitalPhi]d\[Phi]},
    \[Omega] = m \[CapitalOmega];
    rPlus = 1 + Sqrt[1 - a^2];
    iMax = Length[rStarList];
    jMax = Length[thetaList];
    rStar0 = (rStarList[[(iSourceMax + iSourceMin - 1) / 2]] + rStarList[[
      (iSourceMax + iSourceMin + 1) / 2]]) / 2;
    rSubrPlus0 = (rSubrPlusList[[(iSourceMax + iSourceMin - 1) / 2]] 
      + rSubrPlusList[[(iSourceMax + iSourceMin + 1) / 2]]) / 2;
    r0 = rSubrPlus0 + rPlus;
    \[CapitalDelta]rStar = (rStarList[[-1]] - rStarList[[1]]) / (iMax - 1);
    \[CapitalDelta]\[Theta] = (thetaList[[-1]] - thetaList[[1]]) / (jMax - 1);
    solSource = Partition[sol, jMax];
    interp[rStar_, \[Theta]_] = Interpolation[Flatten[Table[{rStarList[[i]],
       thetaList[[j]], solSource[[i, j]]}, {i, (iSourceMax + iSourceMin - 1
      ) / 2 - 1, (iSourceMax + iSourceMin + 1) / 2 + 1}, {j, (jSourceMax + 
      jSourceMin - 1) / 2 - 1, (jSourceMax + jSourceMin + 1) / 2 + 1}], 1],
       InterpolationOrder -> 3][rStar, \[Theta]];
    interpDrStar[rStar_, \[Theta]_] = D[interp[rStar, \[Theta]], rStar];
    \[CapitalPsi] = interp[rStar0, \[Pi] / 2];
    d\[CapitalPsi]drStar = interpDrStar[rStar0, \[Pi] / 2];
    d\[CapitalPhi]dr =
      If[m == 0,
        Re[drStardr[a, rSubrPlus0] d\[CapitalPsi]drStar / r0 - \[CapitalPsi] / r0^2]
        ,
        2 Re[(drStardr[a, rSubrPlus0] d\[CapitalPsi]drStar / r0 + \[CapitalPsi] (I m d\[CapitalDelta]\[Phi]dr[a,
           rSubrPlus0] / r0 - 1 / r0^2)) Exp[I m \[CapitalDelta]\[Phi][a, rSubrPlus0]]]
      ];
    d\[CapitalPhi]dt =
      If[m == 0,
        0
        ,
        2 Re[-I \[Omega] \[CapitalPsi] Exp[I m \[CapitalDelta]\[Phi][a, rSubrPlus0]]] / r0
      ];
    d\[CapitalPhi]d\[Phi] =
      If[m == 0,
        0
        ,
        2 Re[I m \[CapitalPsi] Exp[I m \[CapitalDelta]\[Phi][a, rSubrPlus0]]] / r0
      ];
    {d\[CapitalPhi]dt, d\[CapitalPhi]dr, d\[CapitalPhi]d\[Phi]}
  ]


mRunField[m_, wtDiam_, thetaSourceSize_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, rStarHguess_, 
  rStarIguess_] :=
  Module[{sol, rStarList, rSubrPlusList, thetaList, iSourceMin, iSourceMax,
     jSourceMin, jSourceMax, iMax, jMax, mat, vec, rStarH, rStarL, rStar0,
     rStarR, rStarI},
    {rStarH, rStarL, rStar0, rStarR, rStarI} = getrStarParams[a, r0, 
      wtDiam, rStarHguess, rStarIguess];
    rStarList = getrStarList[\[CapitalDelta]rStar, rStarH, rStarL, rStar0, rStarR, 
      rStarI];
    rSubrPlusList = getrSubrPlusList[a, rStarList];
    thetaList = getThetaList[\[CapitalDelta]\[Theta], thetaSourceSize];
    {iSourceMin, iSourceMax} = getIsourceBounds[rStarList, rStarH, rStarL,
       rStar0, rStarR, rStarI];
    {jSourceMin, jSourceMax} = getJsourceBounds[\[CapitalDelta]\[Theta], thetaSourceSize];
      
    iMax = Length[rStarList];
    jMax = Length[thetaList];
    mat = couplingMatrixFast[m, rStarList, rSubrPlusList, thetaList];
      
    vec = sourceVector[m, rStarList, rSubrPlusList, thetaList, iSourceMin,
       iSourceMax, jSourceMin, jSourceMax];
    sol = LinearSolve[mat, vec, Method -> "Pardiso"];
    sol = Partition[sol, jMax];
    sol = Flatten[Table[{rStarList[[i]],thetaList[[j]],sol[[i,j]]},{i,1,iMax},{j,1,jMax}],1]
  ]


mRunFieldRet[m_, wtDiam_, thetaSourceSize_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, rStarHguess_, 
  rStarIguess_] :=
  Module[{sol, sing, rPlus, rStarList, rSubrPlusList, thetaList, iSourceMin, iSourceMax,
     jSourceMin, jSourceMax, iMax, jMax, mat, vec, rStarH, rStarL, rStar0,
     rStarR, rStarI},
    {rStarH, rStarL, rStar0, rStarR, rStarI} = getrStarParams[a, r0, 
      wtDiam, rStarHguess, rStarIguess];
    rStarList = getrStarList[\[CapitalDelta]rStar, rStarH, rStarL, rStar0, rStarR, 
      rStarI];
    rSubrPlusList = getrSubrPlusList[a, rStarList];
    thetaList = getThetaList[\[CapitalDelta]\[Theta], thetaSourceSize];
    {iSourceMin, iSourceMax} = getIsourceBounds[rStarList, rStarH, rStarL,
       rStar0, rStarR, rStarI];
    {jSourceMin, jSourceMax} = getJsourceBounds[\[CapitalDelta]\[Theta], thetaSourceSize];
      
    iMax = Length[rStarList];
    jMax = Length[thetaList];
    mat = couplingMatrixFast[m, rStarList, rSubrPlusList, thetaList];
      
    vec = sourceVector[m, rStarList, rSubrPlusList, thetaList, iSourceMin,
       iSourceMax, jSourceMin, jSourceMax];
    sol = LinearSolve[mat, vec, Method -> "Pardiso"];
    sol = Partition[sol, jMax];
    rPlus = 1+Sqrt[1-a^2];
    sing = Table[If[(iSourceMin<=i<=iSourceMax)&&(jSourceMin<=j<=jSourceMax),
    PsiPm[m, rSubrPlusList[[i]]+rPlus, thetaList[[j]]]Exp[-I m \[CapitalDelta]\[Phi][a, rSubrPlusList[[i]]]]/(2 \[Pi])
    ,0],{i,1,iMax},{j,1,jMax}];
    sol = Flatten[Table[{rStarList[[i]],thetaList[[j]],sol[[i,j]]+sing[[i,j]]},{i,1,iMax},{j,1,jMax}],1]
  ]


mRun[m_, wtDiam_, thetaSourceSize_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, rStarHguess_, rStarIguess_
  ] :=
  Module[{sol, rStarList, rSubrPlusList, thetaList, iSourceMin, iSourceMax,
     jSourceMin, jSourceMax, iMax, jMax, mat, vec, rStarH, rStarL, rStar0,
     rStarR, rStarI},
    {rStarH, rStarL, rStar0, rStarR, rStarI} = getrStarParams[a, r0, 
      wtDiam, rStarHguess, rStarIguess];
    rStarList = getrStarList[\[CapitalDelta]rStar, rStarH, rStarL, rStar0, rStarR, 
      rStarI];
    rSubrPlusList = getrSubrPlusList[a, rStarList];
    thetaList = getThetaList[\[CapitalDelta]\[Theta], thetaSourceSize];
    {iSourceMin, iSourceMax} = getIsourceBounds[rStarList, rStarH, rStarL,
       rStar0, rStarR, rStarI];
    {jSourceMin, jSourceMax} = getJsourceBounds[\[CapitalDelta]\[Theta], thetaSourceSize];
      
    iMax = Length[rStarList];
    jMax = Length[thetaList];
    mat = couplingMatrixFast[m, rStarList, rSubrPlusList, thetaList];
      
    vec = sourceVector[m, rStarList, rSubrPlusList, thetaList, iSourceMin,
       iSourceMax, jSourceMin, jSourceMax];
    sol = LinearSolve[mat, vec, Method -> "Pardiso"];
    getFm[sol, m, rStarList, rSubrPlusList, thetaList, iSourceMin, iSourceMax,
       jSourceMin, jSourceMax]
  ]


mRunFourth[m_, wtDiam_, thetaSourceSize_, \[CapitalDelta]rStar_, \[CapitalDelta]\[Theta]_, rStarHguess_, rStarIguess_
  ] :=
  Module[{sol, rStarList, rSubrPlusList, thetaList, iSourceMin, iSourceMax,
     jSourceMin, jSourceMax, iMax, jMax, mat, vec, rStarH, rStarL, rStar0,
     rStarR, rStarI},
    {rStarH, rStarL, rStar0, rStarR, rStarI} = getrStarParams[a, r0, 
      wtDiam, rStarHguess, rStarIguess];
    rStarList = getrStarList[\[CapitalDelta]rStar, rStarH, rStarL, rStar0, rStarR, 
      rStarI];
    rSubrPlusList = getrSubrPlusList[a, rStarList];
    thetaList = getThetaList[\[CapitalDelta]\[Theta], thetaSourceSize];
    {iSourceMin, iSourceMax} = getIsourceBounds[rStarList, rStarH, rStarL,
       rStar0, rStarR, rStarI];
    {jSourceMin, jSourceMax} = getJsourceBounds[\[CapitalDelta]\[Theta], thetaSourceSize];
      
    iMax = Length[rStarList];
    jMax = Length[thetaList];
    mat = couplingMatrixFourth[m, rStarList, rSubrPlusList, thetaList];
      
    vec = sourceVectorFourth[m, rStarList, rSubrPlusList, thetaList, iSourceMin,
       iSourceMax, jSourceMin, jSourceMax];
    sol = LinearSolve[mat, vec, Method -> "Pardiso"];
    getFm[sol, m, rStarList, rSubrPlusList, thetaList, iSourceMin, iSourceMax,
       jSourceMin, jSourceMax]
  ]


R2[n_, h1_, h2_, A1_, A2_] := (A2 h1^n - A1 h2^n)/(h1^n - h2^n)


R3[n_, h1_, h2_, h3_, A1_, A2_, A3_] := (A3*h1^n*(h1 - h2)*h2^n + (-(A2*h1^n*(h1 - h3)) + A1*h2^n*(h2 - h3))*h3^n)/
 (h2^n*(h2 - h3)*h3^n + h1^(1 + n)*(h2^n - h3^n) + h1^n*(-h2^(1 + n) + h3^(1 + n)))


R4[n_, h1_, h2_, h3_, h4_, A1_, A2_, A3_, A4_] := (A4*h1^n*(h1 - h2)*h2^n*(h1 - h3)*(h2 - h3)*h3^n + 
  (-(A3*h1^n*(h1 - h2)*h2^n*(h1 - h4)*(h2 - h4)) + 
    h3^n*(A2*h1^n*(h1 - h3)*(h1 - h4) - A1*h2^n*(h2 - h3)*(h2 - h4))*(h3 - h4))*
   h4^n)/(-(h2^n*(h2 - h3)*h3^n*(h2 - h4)*(h3 - h4)*h4^n) + 
  h1^(2 + n)*(h3^n*(h3 - h4)*h4^n + h2^(1 + n)*(h3^n - h4^n) + 
    h2^n*(-h3^(1 + n) + h4^(1 + n))) + 
  h1^(1 + n)*(h3^n*h4^n*(-h3^2 + h4^2) + h2^(2 + n)*(-h3^n + h4^n) + 
    h2^n*(h3^(2 + n) - h4^(2 + n))) + h1^n*(h3^(1 + n)*(h3 - h4)*h4^(1 + n) + 
    h2^(2 + n)*(h3^(1 + n) - h4^(1 + n)) + h2^(1 + n)*(-h3^(2 + n) + h4^(2 + n))))


resolutionConvergeFm[m_, rStarHguess_, rStarIguess_, wtDiam_, thetaSourceSize_,
   tol_, cumulF_] :=
  Module[{converge, converge2, converge3, \[CapitalDelta]rStarGuess, rStarH, rStarL,
     rStar0, rStarR, rStarI, j, \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], Fm, test, newFm, oldFm},
    {rStarH, rStarL, rStar0, rStarR, rStarI} = getrStarParams[a, r0, 
      wtDiam, rStarHguess, rStarIguess];
    converge = {};
    \[CapitalDelta]rStarGuess = 0.7;
    j = Round[wtDiam / \[CapitalDelta]rStarGuess];
    If[EvenQ[j],
      j = j + 1
    ];
    \[CapitalDelta]rStar = wtDiam / j;
    \[CapitalDelta]\[Theta] = thetaSourceSize / (j);
    Fm = mRun[m, wtDiam, thetaSourceSize, \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], rStarHguess, rStarIguess
      ];
    converge = Append[converge, {\[CapitalDelta]rStar, Fm}];
    \[CapitalDelta]rStarGuess = \[CapitalDelta]rStarGuess * 0.7;
    j = Round[wtDiam / \[CapitalDelta]rStarGuess];
    If[EvenQ[j],
      j = j + 1
    ];
    \[CapitalDelta]rStar = wtDiam / j;
    \[CapitalDelta]\[Theta] = thetaSourceSize / (j);
    Fm = mRun[m, wtDiam, thetaSourceSize, \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], rStarHguess, rStarIguess
      ];
    converge = Append[converge, {\[CapitalDelta]rStar, Fm}];
    oldFm = Fm;
    newFm = R2[2, converge[[-2, 1]], converge[[-1, 1]], converge[[-2,
       2]], converge[[-1, 2]]];
    test = 0;
    While[
      test == 0 ||
        AnyTrue[
          Table[
            Abs[newFm[[j]] - oldFm[[j]]] >
              0.3 tol
                If[Max[Abs[cumulF[[j]]], Abs[newFm[[j]]]] < 10 ^ -10,
                  
                  10^10
                  ,
                  Max[Abs[cumulF[[j]]], Abs[newFm[[j]]]]
                ]
            ,
            {j, 2, 2}
          ]
          ,
          TrueQ
        ]
      ,
      test = 1;
      \[CapitalDelta]rStarGuess = \[CapitalDelta]rStarGuess * 0.7;
      j = Round[wtDiam / \[CapitalDelta]rStarGuess];
      If[EvenQ[j],
        j = j + 1
      ];
      \[CapitalDelta]rStar = wtDiam / j;
      \[CapitalDelta]\[Theta] = thetaSourceSize / (j);
      Fm = mRun[m, wtDiam, thetaSourceSize, \[CapitalDelta]rStar, \[CapitalDelta]\[Theta], rStarHguess, 
        rStarIguess];
      converge = Append[converge, {\[CapitalDelta]rStar, Fm}];
      oldFm = newFm;
      newFm =  R3[2, converge[[-3, 1]], converge[[-2, 1]], converge[[-1,
       1]], converge[[-3, 2]], converge[[-2, 2]], converge[[-1, 2]]];
      Print["\[CapitalDelta]r* = ", \[CapitalDelta]rStar, ", \[CapitalDelta]\[Theta] = ", \[CapitalDelta]\[Theta], ", Fm = ", NumberForm[newFm,
         10]];
      
    ];
    Return[newFm]
  ]


domainConvergeFm[m_, rStarHguess_, wtDiamGuess_, thetaSourceSize_, tol_,
   cumulF_] :=
  Module[{j, wtDiam, rStar0, converge, \[CapitalDelta]rStarIguess, \[CapitalDelta]rStarI, test, newFm,
     oldFm, rStarIguess, Fm},
    j = Ceiling[(2 \[Pi] / \[CapitalOmega]) / wtDiamGuess];
    wtDiam = (2 \[Pi] / \[CapitalOmega]) / j;
    rStar0 = getrStarFromrSubrPlus[a, r0 - (1 - Sqrt[1 - a^2])];
    rStarIguess = Round[80.0 / wtDiam] wtDiam + wtDiam / 2 + rStar0;
      
    converge = {};
    \[CapitalDelta]rStarIguess = 50.0;
    j = Ceiling[\[CapitalDelta]rStarIguess / (2 \[Pi] / \[CapitalOmega])];
    \[CapitalDelta]rStarI = j (2 \[Pi] / \[CapitalOmega]);
    Fm = resolutionConvergeFm[m, rStarHguess, rStarIguess, wtDiam, thetaSourceSize,
       tol, cumulF];
    converge = Append[converge, {1 / rStarIguess, Fm}];
    Print["r*I = ", rStarIguess, ", Fm = ", NumberForm[Fm, 10]];
    rStarIguess = rStarIguess + \[CapitalDelta]rStarI;
    Fm = resolutionConvergeFm[m, rStarHguess, rStarIguess, wtDiam, thetaSourceSize,
       tol, cumulF];
    converge = Append[converge, {1 / rStarIguess, Fm}];
    oldFm = Fm;
    newFm = R2[3, converge[[-2, 1]], converge[[-1, 1]], converge[[-2,
       2]], converge[[-1, 2]]];
    Print["r*I = ", rStarIguess, ", Fm = ", NumberForm[Fm, 10]];
    rStarIguess = rStarIguess + \[CapitalDelta]rStarI;
    test = 0;
    While[
      test == 0 ||
        AnyTrue[
          Table[
            Abs[newFm[[j]] - oldFm[[j]]] >
              tol
              If[Max[Abs[cumulF[[j]]], Abs[newFm[[j]]]] < 10 ^ -10,
                  10^10
                  ,
                  Max[Abs[cumulF[[j]]], Abs[newFm[[j]]]]
                ]
            ,
            {j, 2, 2}
          ]
          ,
          TrueQ
        ]
      ,
      test = 1;
      Fm = resolutionConvergeFm[m, rStarHguess, rStarIguess, wtDiam, 
        thetaSourceSize, tol, cumulF];
      converge = Append[converge, {1 / rStarIguess, Fm}];
      oldFm = newFm;
      newFm = 
    newFm = R3[3, converge[[-3, 1]], converge[[-2, 1]], converge[[-1,
        1]], converge[[-3, 2]], converge[[-2, 2]], converge[[-1, 2]]];
      Print["r*I = ", rStarIguess, ", Fm = ", NumberForm[newFm, 10]];
        
      rStarIguess = rStarIguess + \[CapitalDelta]rStarI;
      
    ];
    Return[newFm]
  ]


getF[rStarHguess_, wtDiamGuess_, thetaSourceSize_, tol_] :=
  Module[{m, cumulF, FmList, Fm},
    m = 0;
    Print["m = ", m];
    Fm = domainConvergeFm[m, rStarHguess, wtDiamGuess, thetaSourceSize,
       0.1 tol, {0, 0, 0}];
    FmList = {Fm};
    cumulF = Fm;
    m = 1;
    Print["m = ", m];
    Fm = domainConvergeFm[m, rStarHguess, wtDiamGuess, thetaSourceSize,
       0.3 tol, {0, 0, 0}];
    AppendTo[FmList, Fm];
    cumulF = cumulF + Fm;
    m = 2;
    Print["m = ", m];
    Fm = domainConvergeFm[m, rStarHguess, wtDiamGuess, thetaSourceSize,
       1.0 tol, {0, 0, 0}];
    AppendTo[FmList, Fm];
    cumulF = cumulF + Fm;
    Do[
      Print["m = ", m];
      Fm = domainConvergeFm[m, rStarHguess, wtDiamGuess, thetaSourceSize,
         3.0 tol, cumulF];
      AppendTo[FmList, Fm];
      cumulF = cumulF + Fm;
      ,
      {m, 3, 20}];
    Return[FmList]
  ]; 
