(* ::Package:: *)

BeginPackage["EllipticSSF`"]


Delta[a_, rSubRPlus_] :=
  Module[{r, rPlus, rMinus},
    rPlus = 1 + Sqrt[1 - a^2];
    rMinus = 1 - Sqrt[1 - a^2];
    r = rSubRPlus + rPlus;
    rSubRPlus (r - rMinus)
  ]


Sigma[a_, rSubRPlus_, \[Theta]_] :=
  Module[{r, rPlus, rMinus},
    rPlus = 1 + Sqrt[1 - a^2];
    rMinus = 1 - Sqrt[1 - a^2];
    r = rSubRPlus + rPlus;
    Sqrt[(r^2 + a^2)^2 - a^2 Delta[a, rSubRPlus] Sin[\[Theta]]^2]
  ]


\[CapitalDelta]\[Phi][a_, rSubRPlus_] :=
  Module[{rPlus, rMinus, r},
    rPlus = 1 + Sqrt[1 - a^2];
    rMinus = 1 - Sqrt[1 - a^2];
    r = rSubRPlus + rPlus;
    a / (rPlus - rMinus) Log[rSubRPlus / (r - rMinus)]
  ]; 

(* Define dummy functions, will need to re-define in the notebook where package is imported *)


scaled\[CapitalPhi]Pm[m_, r_, \[Theta]_] :=
  None; 

scaledSeffm[m_, r_, \[Theta]_] :=
  None; 

(* Here is how these functions should be re-defined in the notebook where package is imported: *)

(*
scaled\[CapitalPhi]Pm[m_,r_,\[Theta]_]:=Module[{rPlus},rPlus=1+Sqrt[1-a^2];PsiPm[m,r,\[Theta]]Exp[-I m \[CapitalDelta]\[Phi][a,r-rPlus]]/(2 \[Pi] r)];
scaledSeffm[m_,r_,\[Theta]_]:=Module[{rPlus},rPlus=1+Sqrt[1-a^2];Seffm[m,r,\[Theta]]Exp[-I m \[CapitalDelta]\[Phi][a,r-rPlus]]/(2 \[Pi])];
*)


getfxy[r_, \[Theta]_, A_, \[CapitalDelta]rstar_, \[CapitalDelta]\[Theta]_, m_, a_, \[CapitalOmega]_] := (* re-implemented for package friendliness *)Module[
  {M},
    M = 1;
    Return[
      If[A == 1,
        (a^2 m^2 r^2 (a^2 + r (-2 M + r)) \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 \[CapitalOmega]^2 + (2 r^2 
          (a^2 + r (-2 M + r)) \[CapitalDelta]rstar^2 + 2 (r^2 (a^2 + r^2)^2 - (a^2 - I a m r
           - M r) (a^2 + r (-2 M + r)) \[CapitalDelta]rstar^2) \[CapitalDelta]\[Theta]^2 + 4 a m^2 M r^3 \[CapitalDelta]rstar^2 
          \[CapitalDelta]\[Theta]^2 \[CapitalOmega] - m^2 r^2 (a^2 + r^2)^2 \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 \[CapitalOmega]^2) Csc[\[Theta]]^2 + m^2 r^2 
          (a^2 + r (-2 M + r)) \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 Csc[\[Theta]]^4) / (r^2 \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 (-a
          ^2 (a^2 + r (-2 M + r)) + (a^2 + r^2)^2 Csc[\[Theta]]^2))
        ,
        If[A == 2,
          (-r (a^2 + r^2)^2 + a (a^3 - I a^2 m r - I m r^3 + a r (-2 
            M + r)) \[CapitalDelta]rstar) / (r \[CapitalDelta]rstar^2 ((a^2 + r^2)^2 - a^2 (a^2 + r (-2 M + r
            )) Sin[\[Theta]]^2))
          ,
          If[A == 3,
            -(r (a^2 + r^2)^2 + a (a^3 - I a^2 m r - I m r^3 + a r (-
              2 M + r)) \[CapitalDelta]rstar) / (r \[CapitalDelta]rstar^2 ((a^2 + r^2)^2 - a^2 (a^2 + r (-2 M +
               r)) Sin[\[Theta]]^2))
            ,
            If[A == 4,
              (a^2 + r (-2 M + r)) (2 + \[CapitalDelta]\[Theta] Cot[\[Theta]]) / (2 \[CapitalDelta]\[Theta]^2 (-(a^2 +
                 r^2)^2 + a^2 (a^2 + r (-2 M + r)) Sin[\[Theta]]^2))
              ,
              -(a^2 + r (-2 M + r)) (-2 + \[CapitalDelta]\[Theta] Cot[\[Theta]]) / (2 \[CapitalDelta]\[Theta]^2 (-(a^2
                 + r^2)^2 + a^2 (a^2 + r (-2 M + r)) Sin[\[Theta]]^2))
            ]
          ]
        ]
      ]
    ]
  ]

(* original version below *)

(*
getfxy[r_,\[Theta]_,A_,\[CapitalDelta]rstar_,\[CapitalDelta]\[Theta]_,m_]:=Module[{\[CapitalDelta],\[CapitalSigma]sq,box,\[Delta]t2,\[Delta]t,\[Delta]rstar,\[Delta]rstar2,\[Delta]\[Theta],\[Delta]\[Theta]2,none,fxy,fxyp,fxym,fxpy,fxmy,\[Omega],returnvalue},\[Omega]=m*\[CapitalOmega];\[CapitalDelta][r,\[Theta]]= r^2-2 M r +a^2;\[CapitalSigma]sq[r,\[Theta]]= (r^2+a^2)^2-a^2 ( r^2-2 M r +a^2) Sin[\[Theta]]^2;box[r,\[Theta]] = (\[Delta]t2+((4 I a m M r)/\[CapitalSigma]sq[r,\[Theta]])\[Delta]t-((r^2+a^2)^2/\[CapitalSigma]sq[r,\[Theta]])\[Delta]rstar2-((2 I a m r (r^2+a^2)-2a^2 \[CapitalDelta][r,\[Theta]])/(r \[CapitalSigma]sq[r,\[Theta]]))\[Delta]rstar-(\[CapitalDelta][r,\[Theta]]/\[CapitalSigma]sq[r,\[Theta]])\[Delta]\[Theta]2-(\[CapitalDelta][r,\[Theta]]/\[CapitalSigma]sq[r,\[Theta]] Cot[\[Theta]])\[Delta]\[Theta]-(\[CapitalDelta][r,\[Theta]]/\[CapitalSigma]sq[r,\[Theta]] (-m^2/Sin[\[Theta]]^2-(2M)/r (1-a^2/(M r))-(2 I a m )/r))none) ;(*getting rid of e^(-\[ImaginaryI] \[Omega] t)*)
\[Delta]t2=-\[Omega]^2 Subscript[\[Psi], ij]; \[Delta]t=-I \[Omega] Subscript[\[Psi], ij]; \[Delta]rstar2=(Subscript[\[Psi], ipj]-2Subscript[\[Psi], ij]+Subscript[\[Psi], imj])/\[CapitalDelta]rstar^2; \[Delta]rstar=(Subscript[\[Psi], ipj]-Subscript[\[Psi], imj])/(2 \[CapitalDelta]rstar); \[Delta]\[Theta]2=(Subscript[\[Psi], ijp]-2Subscript[\[Psi], ij]+Subscript[\[Psi], ijm])/\[CapitalDelta]\[Theta]^2; \[Delta]\[Theta]=(Subscript[\[Psi], ijp]-Subscript[\[Psi], ijm])/(2 \[CapitalDelta]\[Theta]);none = Subscript[\[Psi], ij]; 
fxy[r,\[Theta]]=Coefficient[box[r,\[Theta]],Subscript[\[Psi], ij]]//Simplify;fxpy[r,\[Theta]]=Coefficient[box[r,\[Theta]],Subscript[\[Psi], ipj]]//Simplify;
fxmy[r,\[Theta]]=Coefficient[box[r,\[Theta]],Subscript[\[Psi], imj]]//Simplify;fxyp[r,\[Theta]]=Coefficient[box[r,\[Theta]],Subscript[\[Psi], ijp]]//Simplify;
fxym[r,\[Theta]]=Coefficient[box[r,\[Theta]],Subscript[\[Psi], ijm]]//Simplify;returnvalue =If[A==1,fxy[r,\[Theta]],If[A==2,fxpy[r,\[Theta]],If[A==3,fxmy[r,\[Theta]],If[A==4,fxyp[r,\[Theta]],fxym[r,\[Theta]]]]]];returnvalue]
*)


SourceCombined[n_, m_, a_, r0_, rminapp_, rmaxapp_] := (* re-implemented for package friendliness *)Module[
  {rplus, rminus, rtorstar, r0star, rsmin, rsmax, rstartor, listr, listrsource,
   \[Theta]smin, \[Theta]smax, imax, jmax, rsourcestart, rsourceend, listSeff1, SourceMatrix1,
   i, j, r, LeftBC, RightBC, TopBC, BottomBC, SecondTopBC, SecondBottomBC,
   SecondLeftBC, SourceBC1, SecondRightBC, Sourcelist, \[CapitalDelta]\[Theta], \[Theta]sourcegrid,
   l, \[CapitalDelta]rstar, d, rmatrixfactorminus, rmatrixfactorplus, rsSourceMin, rsSourceMax, \[Theta]start,
   \[Theta]end, M, v, \[CapitalOmega], guess},
    M = 1;
    v = 1 / Sqrt[r0];
    \[CapitalOmega] = v^3 / (1 + a v^3);
    \[CapitalDelta]\[Theta] = \[Pi] / (5 (2 n + 1));
    \[Theta]sourcegrid = 2 n + 1;
    l = 3 n + 1;
    \[CapitalDelta]rstar = d / (2 l + 1);
    d = 6 / 5 3 M;
    rplus = M + Sqrt[M^2 - a^2];
    rminus = M - Sqrt[M^2 - a^2];
    rtorstar[r_] := r + (2 M) / (rplus - rminus) (rplus Log[(r - rplus
      ) / (2 M)] - rminus Log[(r - rminus) / (2 M)]);
    r0star = rtorstar[r0];
    rplus = M + Sqrt[M^2 - a^2];
    rminus = M - Sqrt[M^2 - a^2];
    rtorstar[r_] = r + (2 M) / (rplus - rminus) (rplus Log[(r - rplus
      ) / (2 M)] - rminus Log[(r - rminus) / (2 M)]);
    guess[rstar_] :=
      If[rstar <= -2,
        rplus + 2 * (1 - a^2)^(1 / rplus - 1 / 2) Exp[(a^2 - rplus + 
          (rplus - 1) * rstar) / rplus]
        ,
        If[-2 < rstar <= 1000,
          rplus + 2 (ProductLog[Exp[1 / 2 (rstar - rplus)]])
          ,
          rstar + rplus
        ]
      ];
    rstartor[rstar_] :=
      If[rstar < -50,
        guess[rstar]
        ,
        r /. FindRoot[rtorstar[r] - rstar == 0, {r, guess[rstar]}, Method
           -> "Newton", AccuracyGoal -> 13, PrecisionGoal -> 13]
      ];
    r0star = rtorstar[r0];
    \[Theta]smin = \[Pi] / 2 - \[CapitalDelta]\[Theta] / 2 - n * \[CapitalDelta]\[Theta];
    \[Theta]smax = \[Pi] / 2 + \[CapitalDelta]\[Theta] / 2 + n * \[CapitalDelta]\[Theta];
    rsSourceMin = r0star - \[CapitalDelta]rstar / 2 - l * \[CapitalDelta]rstar;
    rsSourceMax = r0star + \[CapitalDelta]rstar / 2 + l * \[CapitalDelta]rstar;
    rmatrixfactorminus = Floor[(rsSourceMin - rminapp) / d] + 1;
    rmatrixfactorplus = Floor[(rmaxapp - rsSourceMax) / d] + 1;
    rsmin = rsSourceMin - (rmatrixfactorminus - 1) * d;
    rsmax = (rmatrixfactorplus - 1) * d + rsSourceMax;
    imax = Floor[(rsmax - rsmin) / \[CapitalDelta]rstar + 1];
    \[Theta]start = \[Theta]smin / \[CapitalDelta]\[Theta];
    \[Theta]end = \[Theta]smax / \[CapitalDelta]\[Theta];
    jmax = Floor[\[Pi] / \[CapitalDelta]\[Theta] + 1];
    listr = Table[Re[rstartor[rstar]], {rstar, rsmin, rsmax, \[CapitalDelta]rstar}];
    listrsource = Table[Re[rstartor[rstar]], {rstar, rsSourceMin, rsSourceMax, \[CapitalDelta]rstar
      }];
    rsourcestart = Floor[(rsSourceMin - rsmin) / \[CapitalDelta]rstar + 1];
    rsourceend = rsourcestart + Length[listrsource] - 1;
    listSeff1 = Table[scaledSeffm[m, listrsource[[i]], N[j]], {i, 1, 
      Length[listrsource]}, {j, \[Theta]smin, \[Theta]smax, \[CapitalDelta]\[Theta]}];
    SourceMatrix1 = SparseArray[Flatten[Table[{i + rsourcestart - 1, 
      j + \[Theta]start} -> listSeff1[[i, j]], {i, 1, Length[listrsource]}, {j, 1,
       \[Theta]end - \[Theta]start + 1}]], {imax, jmax}];
    LeftBC := SparseArray[Flatten[Table[{i2 + rsourcestart - 1, \[Theta]start
      } -> -getfxy[listrsource[[i2]], \[Theta]smin - \[CapitalDelta]\[Theta], 4, \[CapitalDelta]rstar, \[CapitalDelta]\[Theta], m, a, \[CapitalOmega]] listrsource[[
      i2]] scaled\[CapitalPhi]Pm[m, listrsource[[i2]], N[\[Theta]smin]], {i2, 1, Length[listrsource
      ]}]], {imax, jmax}];
    SecondLeftBC := SparseArray[Flatten[Table[{i2 + rsourcestart - 1,
       \[Theta]start + 1} -> getfxy[listrsource[[i2]], \[Theta]smin, 5, \[CapitalDelta]rstar, \[CapitalDelta]\[Theta], m, a,
       \[CapitalOmega]] listrsource[[i2]] scaled\[CapitalPhi]Pm[m, listrsource[[i2]], N[\[Theta]smin - \[CapitalDelta]\[Theta]]],
       {i2, 1, Length[listrsource]}]], {imax, jmax}];
    RightBC := SparseArray[Flatten[Table[{i1 + rsourcestart - 1, \[Theta]end
       + 2} -> -getfxy[listrsource[[i1]], \[Theta]smax + \[CapitalDelta]\[Theta], 5, \[CapitalDelta]rstar, \[CapitalDelta]\[Theta], m, a, 
      \[CapitalOmega]] listrsource[[i1]] scaled\[CapitalPhi]Pm[m, listrsource[[i1]], N[\[Theta]smax]], {i1, 
      1, Length[listrsource]}]], {imax, jmax}];
    SecondRightBC := SparseArray[Flatten[Table[{i1 + rsourcestart - 1,
       \[Theta]end + 1} -> getfxy[listrsource[[i1]], \[Theta]smax, 4, \[CapitalDelta]rstar, \[CapitalDelta]\[Theta], m, a, \[CapitalOmega]
      ] listrsource[[i1]] scaled\[CapitalPhi]Pm[m, listrsource[[i1]], N[\[Theta]smax + \[CapitalDelta]\[Theta]]], {
      i1, 1, Length[listrsource]}]], {imax, jmax}];
    TopBC := SparseArray[Flatten[Table[{rsourcestart - 1, \[Theta]start + j}
       -> -getfxy[listr[[rsourcestart - 1]], \[Theta]smin + \[CapitalDelta]\[Theta] * (j - 1), 2, \[CapitalDelta]rstar,
       \[CapitalDelta]\[Theta], m, a, \[CapitalOmega]] listr[[rsourcestart]] scaled\[CapitalPhi]Pm[m, listr[[rsourcestart]
      ], N[\[Theta]smin + \[CapitalDelta]\[Theta] * (j - 1)]], {j, 1, \[Theta]end - \[Theta]start + 1}]], {imax, jmax
      }];
    SecondTopBC := SparseArray[Flatten[Table[{rsourcestart, \[Theta]start + 
      j} -> getfxy[listr[[rsourcestart]], \[Theta]smin + \[CapitalDelta]\[Theta] * (j - 1), 3, \[CapitalDelta]rstar, 
      \[CapitalDelta]\[Theta], m, a, \[CapitalOmega]] listr[[rsourcestart - 1]] scaled\[CapitalPhi]Pm[m, listr[[rsourcestart
       - 1]], N[\[Theta]smin + \[CapitalDelta]\[Theta] * (j - 1)]], {j, 1, \[Theta]end - \[Theta]start + 1}]], {imax,
       jmax}];
    BottomBC := SparseArray[Flatten[Table[{rsourceend + 1, \[Theta]start + j
      } -> -getfxy[listr[[rsourceend + 1]], \[Theta]smin + \[CapitalDelta]\[Theta] * (j - 1), 3, \[CapitalDelta]rstar,
       \[CapitalDelta]\[Theta], m, a, \[CapitalOmega]] listr[[rsourceend]] scaled\[CapitalPhi]Pm[m, listr[[rsourceend]], N[
      \[Theta]smin + \[CapitalDelta]\[Theta] * (j - 1)]], {j, 1, \[Theta]end - \[Theta]start + 1}]], {imax, jmax}];
    SecondBottomBC := SparseArray[Flatten[Table[{rsourceend, \[Theta]start +
       j} -> getfxy[listr[[rsourceend]], \[Theta]smin + \[CapitalDelta]\[Theta] * (j - 1), 2, \[CapitalDelta]rstar, \[CapitalDelta]\[Theta],
       m, a, \[CapitalOmega]] listr[[rsourceend + 1]] scaled\[CapitalPhi]Pm[m, listr[[rsourceend + 1]
      ], N[\[Theta]smin + \[CapitalDelta]\[Theta] * (j - 1)]], {j, 1, \[Theta]end - \[Theta]start + 1}]], {imax, jmax
      }];
    SourceBC1 = LeftBC + SourceMatrix1 + SecondLeftBC + RightBC + SecondRightBC
       + TopBC + SecondTopBC + BottomBC + SecondBottomBC;
    Sourcelist = Flatten[ArrayReshape[SourceBC1, {1, imax * jmax}]];
    Return[Sourcelist]
  ]

(* original version: SourceCombined[n_,m_] := ... *)


CouplingMatrix[n_, m_, a_, r0_, rminapp_, rmaxapp_] := Module[
  {rplus, rminus, rtorstar, r0star, rstarmin1, rstarmax1, rsmin1, rsmax1,
   \[Theta]smin, \[Theta]smax, imax, jmax, listr1, listrsource1, MxBC, rstartor, listr,
   listrsource, xBC1, xBC2, xBC3, yBC1, yBC2, yBC3, \[Omega], d, rsmin, rsmax,
   rsSourceMin, rsSourceMax, \[CapitalDelta]\[Theta], \[Theta]sourcegrid, l, \[CapitalDelta]rstar, rsourcegrid, rmatrixfactorminus,
   rmatrixfactorplus, M, v, \[CapitalOmega], listTheta, mat1, mat2, diag, rightDiag, 
  leftDiag, rightSkipDiag, leftSkipDiag, guess},
    M = 1;
    v = 1 / Sqrt[r0];
    \[CapitalOmega] = v^3 / (1 + a v^3);
    \[Omega] = m * \[CapitalOmega];
    \[CapitalDelta]\[Theta] = \[Pi] / (5 (2 n + 1));
    \[Theta]sourcegrid = 2 n + 1;
    l = 3 n + 1;
    \[CapitalDelta]rstar = d / (2 l + 1);
    rsourcegrid = 2 l + 1;
    d = 6 / 5 3 M;
    rplus = M + Sqrt[M^2 - a^2];
    rminus = M - Sqrt[M^2 - a^2];
    rtorstar[r_] := r + (2 M) / (rplus - rminus) (rplus Log[(r - rplus
      ) / (2 M)] - rminus Log[(r - rminus) / (2 M)]);
    r0star = rtorstar[r0];
    \[Theta]smin = \[Pi] / 2 - \[CapitalDelta]\[Theta] / 2 - n * \[CapitalDelta]\[Theta];
    \[Theta]smax = \[Pi] / 2 + \[CapitalDelta]\[Theta] / 2 + n * \[CapitalDelta]\[Theta];
    rsSourceMin = r0star - \[CapitalDelta]rstar / 2 - l * \[CapitalDelta]rstar;
    rsSourceMax = r0star + \[CapitalDelta]rstar / 2 + l * \[CapitalDelta]rstar;
    rmatrixfactorminus = Floor[(rsSourceMin - rminapp) / d] + 1;
    rmatrixfactorplus = Floor[(rmaxapp - rsSourceMax) / d] + 1;
    rsmin = rsSourceMin - (rmatrixfactorminus - 1) * d;
    rsmax = (rmatrixfactorplus - 1) * d + rsSourceMax;
    imax = Floor[(rsmax - rsmin) / \[CapitalDelta]rstar];
    jmax = Floor[\[Pi] / \[CapitalDelta]\[Theta]];
    guess[rstar_] :=
      If[rstar <= -2,
        rplus + 2 * (1 - a^2)^(1 / rplus - 1 / 2) Exp[(a^2 - rplus + 
          (rplus - 1) * rstar) / rplus]
        ,
        If[-2 < rstar <= 1000,
          rplus + 2 (ProductLog[Exp[1 / 2 (rstar - rplus)]])
          ,
          rstar + rplus
        ]
      ];
    rstartor[rstar_] :=
      If[rstar < -50,
        guess[rstar]
        ,
        r /. FindRoot[rtorstar[r] - rstar == 0, {r, guess[rstar]}, Method
           -> "Newton", AccuracyGoal -> 13, PrecisionGoal -> 13]
      ];
    listr1 = Table[Re[rstartor[rstar]], {rstar, rsmin + \[CapitalDelta]rstar,
       rsmax - \[CapitalDelta]rstar, \[CapitalDelta]rstar}];
    listTheta = Table[j * \[CapitalDelta]\[Theta], {j, 1, jmax - 1}];
    xBC1 = 3.0 - 2 I \[CapitalDelta]rstar \[Omega];
    xBC2 = -4.0;
    xBC3 = 1.0;
    yBC1 =
      If[m == 0,
        3.0
        ,
        1.0
      ];
    yBC2 =
      If[m == 0,
        -4.0
        ,
        0
      ];
    yBC3 =
      If[m == 0,
        1.0
        ,
        0
      ];
    diag = Flatten[ArrayPad[Transpose[Table[(a^2 m^2 listr1^2
       (a^2 + listr1 (-2 M + listr1)) \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 \[CapitalOmega]^2 + (2 listr1^2 (a^2 
      + listr1 (-2 M + listr1)) \[CapitalDelta]rstar^2 + 2 (listr1^2 (a^2 + listr1^2)^2 -
       (a^2 - I a m listr1 - M listr1) (a^2 + listr1 (-2 M + listr1)) \[CapitalDelta]rstar
      ^2) \[CapitalDelta]\[Theta]^2 + 4 a m^2 M listr1^3 \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 \[CapitalOmega] - m^2 listr1^2 (a^2 + listr1
      ^2)^2 \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 \[CapitalOmega]^2) Csc[theta]^2 + m^2 listr1^2 (a^2 + listr1 (-
      2 M + listr1)) \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 Csc[theta]^4) / (listr1^2 \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 
      (-a^2 (a^2 + listr1 (-2 M + listr1)) + (a^2 + listr1^2)^2 Csc[theta]^
      2)), {theta, listTheta}]], 1]];
    rightDiag = Flatten[ArrayPad[Transpose[Table[(a^2 + listr1
       (-2 M + listr1)) (2 + \[CapitalDelta]\[Theta] Cot[theta]) / (2 \[CapitalDelta]\[Theta]^2 (-(a^2 + listr1^2)^2 
      + a^2 (a^2 + listr1 (-2 M + listr1)) Sin[theta]^2)), {theta, listTheta
      }]], 1]];
    leftDiag = Flatten[ArrayPad[Transpose[Table[-(a^2 + listr1
       (-2 M + listr1)) (-2 + \[CapitalDelta]\[Theta] Cot[theta]) / (2 \[CapitalDelta]\[Theta]^2 (-(a^2 + listr1^2)^2
       + a^2 (a^2 + listr1 (-2 M + listr1)) Sin[theta]^2)), {theta, listTheta
      }]], 1]];
    rightSkipDiag = Flatten[ArrayPad[Transpose[Table[(-listr1
       (a^2 + listr1^2)^2 + a (a^3 - I a^2 m listr1 - I m listr1^3 + a listr1
       (-2 M + listr1)) \[CapitalDelta]rstar) / (listr1 \[CapitalDelta]rstar^2 ((a^2 + listr1^2)^2 - a^
      2 (a^2 + listr1 (-2 M + listr1)) Sin[theta]^2)), {theta, listTheta}]],
       {{1, 0}, {1, 1}}]];
    leftSkipDiag = Flatten[ArrayPad[Transpose[Table[-(listr1 
      (a^2 + listr1^2)^2 + a (a^3 - I a^2 m listr1 - I m listr1^3 + a listr1
       (-2 M + listr1)) \[CapitalDelta]rstar) / (listr1 \[CapitalDelta]rstar^2 ((a^2 + listr1^2)^2 - a^
      2 (a^2 + listr1 (-2 M + listr1)) Sin[theta]^2)), {theta, listTheta}]],
       {{0, 1}, {1, 1}}]];
    MxBC = SparseArray[{Band[{1, 1}, {(jmax + 1), (jmax + 1)}] -> xBC1,
       Band[{1 + imax + imax jmax, 1 + imax + imax jmax}, {(imax + 1) (jmax
       + 1), (imax + 1) (jmax + 1)}] -> xBC1, Band[{1, jmax + 2}, {(jmax + 
      1), 2 (1 + jmax)}] -> xBC2, Band[{1 + imax + imax jmax, imax - jmax +
       imax jmax}, {(imax + 1) (jmax + 1), imax (1 + jmax)}] -> xBC2, Band[
      {1, 3 + 2 jmax}, {(jmax + 1), 3 (1 + jmax)}] -> xBC3, Band[{1 + imax 
      + imax jmax, -1 + imax - 2 jmax + imax jmax}, {(imax + 1) (jmax + 1),
       (-1 + imax) (1 + jmax)}] -> xBC3, Band[{2 + jmax, 2 + jmax}, {imax +
       (-1 + imax) jmax, imax + (-1 + imax) jmax}, 1 + jmax] -> yBC1, Band[
      {2 (1 + jmax), 2 (1 + jmax)}, {(1 + imax - 1) (1 + jmax), (1 + imax -
       1) (1 + jmax)}, 1 + jmax] -> yBC1, Band[{2 + jmax, 3 + jmax}, {imax 
      + (-1 + imax) jmax, imax + (-1 + imax) jmax + 1}, 1 + jmax] -> yBC2, 
      Band[{2 (1 + jmax), 2 (1 + jmax) - 1}, {(1 + imax - 1) (1 + jmax), (1
       + imax - 1) (1 + jmax) - 1}, 1 + jmax] -> yBC2, Band[{2 + jmax, 4 + 
      jmax}, {imax + (-1 + imax) jmax, imax + (-1 + imax) jmax + 2}, 1 + jmax
      ] -> yBC3, Band[{2 (1 + jmax), 2 (1 + jmax) - 2}, {(1 + imax - 1) (1 
      + jmax), (1 + imax - 1) (1 + jmax) - 2}, 1 + jmax] -> yBC3}, {(imax +
       1) (jmax + 1), (imax + 1) (jmax + 1)}, 0] + DiagonalMatrix[SparseArray[
      {i_} :> diag[[i]], {Length[diag]}, 0], 0] + DiagonalMatrix[SparseArray[
      {i_} :> rightDiag[[i]], {Length[rightDiag] - 1}, 0], 1] + DiagonalMatrix[
      SparseArray[{i_} :> leftDiag[[i + 1]], {Length[leftDiag] - 1}, 0], -1
      ] + DiagonalMatrix[SparseArray[{i_} :> rightSkipDiag[[i]], {Length[rightSkipDiag
      ]}, 0], jmax + 1] + DiagonalMatrix[SparseArray[{i_} :> leftSkipDiag[[
      i]], {Length[leftSkipDiag]}, 0], -jmax - 1];
    Return[MxBC]
  ]


CouplingMatrixOld[n_, m_, a_, r0_, rminapp_, rmaxapp_] := (* re-implemented for package-friendliness *)Module[
  {rplus, rminus, rtorstar, r0star, rstarmin1, rstarmax1, rsmin1, rsmax1,
   \[Theta]smin, \[Theta]smax, imax, jmax, listr1, listrsource1, MxBC, rstartor, listr,
   listrsource, xBC1, xBC2, xBC3, yBC1, yBC2, yBC3, \[Omega], d, rsmin,
   rsmax, rsSourceMin, rsSourceMax, \[CapitalDelta]\[Theta], \[Theta]sourcegrid, l, \[CapitalDelta]rstar, rsourcegrid, rmatrixfactorminus,
   rmatrixfactorplus, M, v, \[CapitalOmega], guess},
    M = 1;
    v = 1 / Sqrt[r0];
    \[CapitalOmega] = v^3 / (1 + a v^3);
    \[Omega] = m * \[CapitalOmega];
    \[CapitalDelta]\[Theta] = \[Pi] / (5 (2 n + 1));
    \[Theta]sourcegrid = 2 n + 1;
    l = 3 n + 1;
    \[CapitalDelta]rstar = d / (2 l + 1);
    rsourcegrid = 2 l + 1;
    d = 6 / 5 3 M;
    rplus = M + Sqrt[M^2 - a^2];
    rminus = M - Sqrt[M^2 - a^2];
    rtorstar[r_] := r + (2 M) / (rplus - rminus) (rplus Log[(r - rplus
      ) / (2 M)] - rminus Log[(r - rminus) / (2 M)]);
    r0star = rtorstar[r0];
    \[Theta]smin = \[Pi] / 2 - \[CapitalDelta]\[Theta] / 2 - n * \[CapitalDelta]\[Theta];
    \[Theta]smax = \[Pi] / 2 + \[CapitalDelta]\[Theta] / 2 + n * \[CapitalDelta]\[Theta];
    rsSourceMin = r0star - \[CapitalDelta]rstar / 2 - l * \[CapitalDelta]rstar;
    rsSourceMax = r0star + \[CapitalDelta]rstar / 2 + l * \[CapitalDelta]rstar;
    rmatrixfactorminus = Floor[(rsSourceMin - rminapp) / d] + 1;
    rmatrixfactorplus = Floor[(rmaxapp - rsSourceMax) / d] + 1;
    rsmin = rsSourceMin - (rmatrixfactorminus - 1) * d;
    rsmax = (rmatrixfactorplus - 1) * d + rsSourceMax;
    imax = Floor[(rsmax - rsmin) / \[CapitalDelta]rstar];
    jmax = Floor[\[Pi] / \[CapitalDelta]\[Theta]];
    guess[rstar_] :=
      If[rstar <= -2,
        rplus + 2 * (1 - a^2)^(1 / rplus - 1 / 2) Exp[(a^2 - rplus + 
          (rplus - 1) * rstar) / rplus]
        ,
        If[-2 < rstar <= 1000,
          rplus + 2 (ProductLog[Exp[1 / 2 (rstar - rplus)]])
          ,
          rstar + rplus
        ]
      ];
    rstartor[rstar_] :=
      If[rstar < -50,
        guess[rstar]
        ,
        r /. FindRoot[rtorstar[r] - rstar == 0, {r, guess[rstar]}, Method
           -> "Newton", AccuracyGoal -> 13, PrecisionGoal -> 13]
      ];
    listr1 = Table[Re[rstartor[rstar]], {rstar, rsmin, rsmax, \[CapitalDelta]rstar}];
    xBC1 = 3 - 2 I \[CapitalDelta]rstar \[Omega];
    xBC2 = -4;
    xBC3 = 1;
    yBC1 =
      If[m == 0,
        3
        ,
        1
      ];
    yBC2 =
      If[m == 0,
        -4
        ,
        0
      ];
    yBC3 =
      If[m == 0,
        1
        ,
        0
      ];
    MxBC = SparseArray[Flatten[Table[{Band[{1, 1}, {(jmax + 1), (jmax
       + 1)}] -> xBC1, Band[{1 + imax + imax jmax, 1 + imax + imax jmax}, {
      (imax + 1) (jmax + 1), (imax + 1) (jmax + 1)}] -> xBC1, Band[{1, jmax
       + 2}, {(jmax + 1), 2 (1 + jmax)}] -> xBC2, Band[{1 + imax + imax jmax,
       imax - jmax + imax jmax}, {(imax + 1) (jmax + 1), imax (1 + jmax)}] 
      -> xBC2, Band[{1, 3 + 2 jmax}, {(jmax + 1), 3 (1 + jmax)}] -> xBC3, Band[
      {1 + imax + imax jmax, -1 + imax - 2 jmax + imax jmax}, {(imax + 1) (
      jmax + 1), (-1 + imax) (1 + jmax)}] -> xBC3, Band[{2 + jmax, 2 + jmax
      }, {imax + (-1 + imax) jmax, imax + (-1 + imax) jmax}, 1 + jmax] -> yBC1,
       Band[{2 (1 + jmax), 2 (1 + jmax)}, {(1 + imax - 1) (1 + jmax), (1 + 
      imax - 1) (1 + jmax)}, 1 + jmax] -> yBC1, Band[{2 + jmax, 3 + jmax}, 
      {imax + (-1 + imax) jmax, imax + (-1 + imax) jmax + 1}, 1 + jmax] -> 
      yBC2, Band[{2 (1 + jmax), 2 (1 + jmax) - 1}, {(1 + imax - 1) (1 + jmax
      ), (1 + imax - 1) (1 + jmax) - 1}, 1 + jmax] -> yBC2, Band[{2 + jmax,
       4 + jmax}, {imax + (-1 + imax) jmax, imax + (-1 + imax) jmax + 2}, 1
       + jmax] -> yBC3, Band[{2 (1 + jmax), 2 (1 + jmax) - 2}, {(1 + imax -
       1) (1 + jmax), (1 + imax - 1) (1 + jmax) - 2}, 1 + jmax] -> yBC3, Band[
      {(3 + jmax) + (jmax + 1) (i - 1), (3 + jmax) + (jmax + 1) (i - 1)}, {
      (2 + jmax) + (jmax + 1) (i - 1) + (jmax - 1), (2 + jmax) + (jmax + 1)
       (i - 1) + (jmax - 1)}] -> Table[getfxy[listr1[[i + 1]], j * \[CapitalDelta]\[Theta], 1, \[CapitalDelta]rstar,
       \[CapitalDelta]\[Theta], m, a, \[CapitalOmega]], {j, 1, jmax - 1}], Band[{(3 + jmax) + (jmax + 1) (i - 
      1), (3 + jmax) + (jmax + 1) (i - 1) - jmax - 1}, {(2 + jmax) + (jmax 
      + 1) (i - 1) + (jmax - 1), (2 + jmax) + (jmax + 1) (i - 1) + (jmax - 
      1) - jmax - 1}] -> Table[getfxy[listr1[[i + 1]], j * \[CapitalDelta]\[Theta], 3, \[CapitalDelta]rstar, \[CapitalDelta]\[Theta],
       m, a, \[CapitalOmega]], {j, 1, jmax - 1}], Band[{(3 + jmax) + (jmax + 1) (i - 1), 
      (3 + jmax) + (jmax + 1) (i - 1) + jmax + 1}, {(2 + jmax) + (jmax + 1)
       (i - 1) + (jmax - 1), (2 + jmax) + (jmax + 1) (i - 1) + (jmax - 1) +
       jmax + 1}] -> Table[getfxy[listr1[[i + 1]], j * \[CapitalDelta]\[Theta], 2, \[CapitalDelta]rstar, \[CapitalDelta]\[Theta], m,
       a, \[CapitalOmega]], {j, 1, jmax - 1}], Band[{(3 + jmax) + (jmax + 1) (i - 1), (3 
      + jmax) + (jmax + 1) (i - 1) - 1}, {(2 + jmax) + (jmax + 1) (i - 1) +
       (jmax - 1), (2 + jmax) + (jmax + 1) (i - 1) + (jmax - 1) - 1}] -> Table[
      getfxy[listr1[[i + 1]], j * \[CapitalDelta]\[Theta], 5, \[CapitalDelta]rstar, \[CapitalDelta]\[Theta], m, a, \[CapitalOmega]], {j, 1, jmax 
      - 1}], Band[{(3 + jmax) + (jmax + 1) (i - 1), (3 + jmax) + (jmax + 1)
       (i - 1) + 1}, {(2 + jmax) + (jmax + 1) (i - 1) + (jmax - 1), (2 + jmax
      ) + (jmax + 1) (i - 1) + (jmax - 1) + 1}] -> Table[getfxy[listr1[[i +
       1]], j * \[CapitalDelta]\[Theta], 4, \[CapitalDelta]rstar, \[CapitalDelta]\[Theta], m, a, \[CapitalOmega]], {j, 1, jmax - 1}]}, {i, 1, imax
       - 1}]], {(imax + 1) (jmax + 1), (imax + 1) (jmax + 1)}, 0];
    Return[MxBC]
  ]

(* original version: CouplingMatrix[n_,m_] := ... *)


\[Psi]FDA[n_, m_, a_, r0_, rminapp_, rmaxapp_] := (* re-implemented for package-friendliness *)Module[
  {final\[Psi]1, jmax, \[CapitalDelta]\[Theta]},
    \[CapitalDelta]\[Theta] = \[Pi] / (5 (2 n + 1));
    jmax = Floor[\[Pi] / \[CapitalDelta]\[Theta] + 1];
    final\[Psi]1 = Partition[LinearSolve[CouplingMatrix[n, m, a, r0, rminapp,
       rmaxapp], SourceCombined[n, m, a, r0, rminapp, rmaxapp], Method -> "Pardiso"
      ], jmax];
    Return[final\[Psi]1]
  ]

(* original version: \[Psi]FDA[n_,m_] := ... *)


\[Psi]atr\[Theta][n_, m_, a_, r0_, rminapp_, rmaxapp_] := (* re-implemented for package-friendliness *)Module[
  {imax, rstarmax, rstarmin, r0star, jmax, rplus, rminus, rtorstar, listr,
   rstartor, listrstar, list\[Theta], function\[CapitalPsi], d\[CapitalPsi]r, \[CurlyPhi], \[CapitalDelta], \[Psi]inr\[Theta]1try1,
   \[CapitalDelta]rstar, \[CapitalDelta]\[Theta], d, rmatrixfactorminus, rmatrixfactorplus, rsmin, rsmax,
  l, M, v, \[CapitalOmega], \[Psi]values, rsSourceMin, rsSourceMax, guess},
    M = 1;
    v = 1 / Sqrt[r0];
    \[CapitalOmega] = v^3 / (1 + a v^3);
    \[CapitalDelta]\[Theta] = \[Pi] / (5 (2 n + 1));
    l = 3 n + 1;
    d = (6 / 5) 3 M;
    \[CapitalDelta]rstar = d / (2 l + 1);
    rplus = M + Sqrt[M^2 - a^2];
    rminus = M - Sqrt[M^2 - a^2];
    rtorstar[r_] := r + (2 M) / (rplus - rminus) (rplus Log[(r - rplus
      ) / (2 M)] - rminus Log[(r - rminus) / (2 M)]);
    r0star = rtorstar[r0];
    rsSourceMin = r0star - \[CapitalDelta]rstar / 2 - l * \[CapitalDelta]rstar;
    rsSourceMax = r0star + \[CapitalDelta]rstar / 2 + l * \[CapitalDelta]rstar;
    rmatrixfactorminus = Floor[(rsSourceMin - rminapp) / d] + 1;
    rmatrixfactorplus = Floor[(rmaxapp - rsSourceMax) / d] + 1;
    rsmin = rsSourceMin - (rmatrixfactorminus - 1) * d;
    rsmax = (rmatrixfactorplus - 1) * d + rsSourceMax;
    guess[rstar_] :=
      If[rstar <= -2,
        rplus + 2 * (1 - a^2)^(1 / rplus - 1 / 2) Exp[(a^2 - rplus + 
          (rplus - 1) * rstar) / rplus]
        ,
        If[-2 < rstar <= 1000,
          rplus + 2 (ProductLog[Exp[1 / 2 (rstar - rplus)]])
          ,
          rstar + rplus
        ]
      ];
    rstartor[rstar_] :=
      If[rstar < -50,
        guess[rstar]
        ,
        r /. FindRoot[rtorstar[r] - rstar == 0, {r, guess[rstar]}, Method
           -> "Newton", AccuracyGoal -> 13, PrecisionGoal -> 13]
      ];
    imax = Floor[(rsmax - rsmin) / \[CapitalDelta]rstar + 1];
    jmax = Floor[\[Pi] / \[CapitalDelta]\[Theta] + 1];
    listr = Table[Re[rstartor[rstar]], {rstar, rsmin, rsmax, \[CapitalDelta]rstar}] ;
    listrstar = Table[rtorstar[listr[[i]]], {i, 1, Length[listr]}];
    list\[Theta] = Table[(j - 1) \[CapitalDelta]\[Theta], {j, 1, jmax}];
    \[CurlyPhi][r_] := a / (rplus - rminus) Log[(r - rplus) / (r - rminus)];
    \[CapitalDelta][r_] := r^2 - 2 M r + a^2;
    \[Psi]values = \[Psi]FDA[n, m, a, r0, rminapp, rmaxapp];
    {\[Psi]inr\[Theta]1try1 = Flatten[Table[{listrstar[[i]], list\[Theta][[j]], 2 * \[Psi]values[[
      i, j]]}, {i, 1, imax}, {j, 1, Length[list\[Theta]]}], 1], imax, jmax, \[CapitalDelta]rstar,
       \[CapitalDelta]\[Theta]}
  ]

(* original version: \[Psi]atr\[Theta][n_,m_] := ... *)


ssf[n_, m_, a_, r0_, rminapp_, rmaxapp_] :=
  Module[{\[Psi]valuesatr\[Theta], function\[Psi], r0star, M, \[CapitalDelta]\[Theta], l, d, \[CapitalDelta]rstar, d\[CapitalPsi]r, rplus,
     rminus, d\[Phi]r1, d\[Phi]\[Theta], d\[Phi]\[Phi]},
    M = 1;
    rplus = M + Sqrt[M^2 - a^2];
    rminus = M - Sqrt[M^2 - a^2];
    \[CapitalDelta]\[Theta] = \[Pi] / (5 (2 n + 1));
    l = 3 n + 1;
    d = (6 / 5) 3 M;
    \[CapitalDelta]rstar = d / (2 l + 1);
    \[Psi]valuesatr\[Theta] = \[Psi]atr\[Theta][n, m, a, r0, rminapp, rmaxapp];
    function\[Psi] = Interpolation[\[Psi]valuesatr\[Theta][[1]]];
    r0star = r0 + (2 M) / (rplus - rminus) (rplus Log[(r0 - rplus) / 
      (2 M)] - rminus Log[(r0 - rminus) / (2 M)]);
    d\[CapitalPsi]r = (-function\[Psi][r0star + 2 \[CapitalDelta]rstar, \[Pi] / 2] + 8 function\[Psi][r0star 
      + \[CapitalDelta]rstar, \[Pi] / 2] - 8 function\[Psi][r0star - \[CapitalDelta]rstar, \[Pi] / 2] + function\[Psi][r0star
       - 2 \[CapitalDelta]rstar, \[Pi] / 2]) / (12 \[CapitalDelta]rstar);
    d\[Phi]r1 = (1 / r0 function\[Psi][r0star, \[Pi] / 2] I m a / ((r0 - rminus) (r0
       - rplus)) + 1 / r0 * (r0^2 + a^2) / Delta[r0] d\[CapitalPsi]r - function\[Psi][r0star,
       \[Pi] / 2] * 1 / r0^2);
    d\[Phi]\[Theta] = (-function\[Psi][r0star, \[Pi] / 2 + 2 \[CapitalDelta]\[Theta]] + 8 function\[Psi][r0star, \[Pi] /
       2 + \[CapitalDelta]\[Theta]] - 8 function\[Psi][r0star, \[Pi] / 2 - \[CapitalDelta]\[Theta]] + function\[Psi][r0star, \[Pi] / 2 
      - 2 \[CapitalDelta]\[Theta]]) / (12 \[CapitalDelta]\[Theta]);
    d\[Phi]\[Phi] = 1 / r0 function\[Psi][r0star, \[Pi] / 2] * (I m);
    Return[{d\[CapitalPsi]r, d\[Phi]r1, d\[Phi]\[Theta], d\[Phi]\[Phi]}]
  ]


EndPackage[]
