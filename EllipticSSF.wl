(* ::Package:: *)

(*Calculating Scalar Self Force Via Elliptic PDEs *)
(*This is a  Mathematica Notebook to calculate the scalar self force via Elliptic PDEs. *)


(*You can change those values to try different mode m, spin parameter a, mass M, and resolutions. *)

m=2; (* mode m *)
n=23; (* resolution *)
rminapp=-20;(*approximate value for the r^* minimum *)
rmaxapp=250;(*approximate value for the r^* maximum  *)
q=1; (* sclar charge *)
d=3.6; (*width of the worldtube *)
a=0.0; (* spin parameter *)
M=1; (* mass of the central black hole *)
r0=10.0;(* location of the worldline *)
v=1/Sqrt[r0];(* particle's angular velocity *)
\[CapitalOmega]=v^3/(1+a v^3);(* particle's angular frequency *)



(*Importing Effective Source  *)

SetDirectory[NotebookDirectory[]]
link=Install["a.out"];
InitializeSeff[r0,a];

(* Here, I'm calling the C-code "kerr-circular.c" created by Barry Wardell, which is available on his website https://github.com/barrywardell/EffectiveSource/blob/master/kerr-circular.c *)
(* "a.out" is this C-code, adjusted to work for Mathematica Notebook *) 
 


(*Effective Source Calculation *)

Delta[a_, rSubRPlus_] :=
  Module[{r, rPlus, rMinus},
    rPlus = 1 + Sqrt[1 - a^2];
    rMinus = 1 - Sqrt[1 - a^2];
    r = rSubRPlus + rPlus;
    rSubRPlus (r - rMinus)
  ];
Sigma[a_, rSubRPlus_, \[Theta]_] :=
  Module[{r, rPlus, rMinus},
    rPlus = 1 + Sqrt[1 - a^2];
    rMinus = 1 - Sqrt[1 - a^2];
    r = rSubRPlus + rPlus;
    Sqrt[(r^2 + a^2)^2 - a^2 Delta[a, rSubRPlus] Sin[\[Theta]]^2]
  ];
\[CapitalDelta]\[Phi][a_, rSubRPlus_] :=
  Module[{rPlus, rMinus, r},
    rPlus = 1 + Sqrt[1 - a^2];
    rMinus = 1 - Sqrt[1 - a^2];
    r = rSubRPlus + rPlus;
    a / (rPlus - rMinus) Log[rSubRPlus / (r - rMinus)]
  ]; 

scaled\[CapitalPhi]Pm[m_,r_,\[Theta]_]:=Module[{rPlus},rPlus=1+Sqrt[1-a^2];PsiPm[m,r,\[Theta]]Exp[-I m \[CapitalDelta]\[Phi][a,r-rPlus]]/(2 \[Pi] r)];
scaledSeffm[m_,r_,\[Theta]_]:=Module[{rPlus},rPlus=1+Sqrt[1-a^2];Seffm[m,r,\[Theta]]Exp[-I m \[CapitalDelta]\[Phi][a,r-rPlus]]/(2 \[Pi])];



(*Solving m-mode scalar wave operator *)

(*Here, I'm solving a m-mode scalar wave operator by using finite difference approximation. 
The coefficient of Subscript[\[Psi], i,j], Subscript[\[Psi], i+1,j], Subscript[\[Psi], i-1, j], Subscript[\[Psi], i,j+1]  and Subscript[\[Psi], i,j-1] are calculated (which are controlled by the parameter "A" inside this function) 
Function below was re-implemented for package friendliness *)

getfxyCoeff[r_, \[Theta]_, A_, \[CapitalDelta]rstar_, \[CapitalDelta]\[Theta]_, m_] :=Module[ {M},
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
  
  
  

(*This is the original function *)

(*getfxyCoeff[r_,\[Theta]_,A_,\[CapitalDelta]rstar_,\[CapitalDelta]\[Theta]_,m_]:=Module[{\[CapitalDelta],\[CapitalSigma]sq,box,\[Delta]t2,\[Delta]t,\[Delta]rstar,\[Delta]rstar2,\[Delta]\[Theta],\[Delta]\[Theta]2,none,fxy,fxyp,fxym,fxpy,fxmy,\[Omega],returnvalue},\[Omega]=m*\[CapitalOmega];\[CapitalDelta][r,\[Theta]]= r^2-2 M r +a^2;\[CapitalSigma]sq[r,\[Theta]]= (r^2+a^2)^2-a^2 ( r^2-2 M r +a^2) Sin[\[Theta]]^2;box[r,\[Theta]] = (\[Delta]t2+((4 \[ImaginaryI] a m M r)/\[CapitalSigma]sq[r,\[Theta]])\[Delta]t-((r^2+a^2)^2/\[CapitalSigma]sq[r,\[Theta]])\[Delta]rstar2-((2 \[ImaginaryI] a m r (r^2+a^2)-2a^2\[CapitalDelta])/(r \[CapitalSigma]sq[r,\[Theta]]))\[Delta]rstar-(\[CapitalDelta][r,\[Theta]]/\[CapitalSigma]sq[r,\[Theta]])\[Delta]\[Theta]2-(\[CapitalDelta][r,\[Theta]]/\[CapitalSigma]sq[r,\[Theta]]Cot[\[Theta]])\[Delta]\[Theta]-(\[CapitalDelta][r,\[Theta]]/\[CapitalSigma]sq[r,\[Theta]](-m^2/Sin[\[Theta]]^2-(2M)/r(1-a^2/(M r))-(2 \[ImaginaryI] a m )/r))none) ;(*getting rid of e^(-\[ImaginaryI] \[Omega] t)*)
\[Delta]t2=-\[Omega]^2Subscript[\[Psi], ij]; \[Delta]t=-\[ImaginaryI] \[Omega] Subscript[\[Psi], ij]; \[Delta]rstar2=(Subscript[\[Psi], ipj]-2Subscript[\[Psi], ij]+Subscript[\[Psi], imj])/\[CapitalDelta]rstar^2; \[Delta]rstar=(Subscript[\[Psi], ipj]-Subscript[\[Psi], imj])/(2 \[CapitalDelta]rstar); \[Delta]\[Theta]2=(Subscript[\[Psi], ijp]-2Subscript[\[Psi], ij]+Subscript[\[Psi], ijm])/\[CapitalDelta]\[Theta]^2; \[Delta]\[Theta]=(Subscript[\[Psi], ijp]-Subscript[\[Psi], ijm])/(2 \[CapitalDelta]\[Theta]);none = Subscript[\[Psi], ij]; 
fxy[r,\[Theta]]=Coefficient[box[r,\[Theta]],Subscript[\[Psi], ij]]//Simplify;fxpy[r,\[Theta]]=Coefficient[box[r,\[Theta]],Subscript[\[Psi], ipj]]//Simplify;
fxmy[r,\[Theta]]=Coefficient[box[r,\[Theta]],Subscript[\[Psi], imj]]//Simplify;fxyp[r,\[Theta]]=Coefficient[box[r,\[Theta]],Subscript[\[Psi], ijp]]//Simplify;
fxym[r,\[Theta]]=Coefficient[box[r,\[Theta]],Subscript[\[Psi], ijm]]//Simplify;returnvalue =If[A\[Equal]1,fxy[r,\[Theta]],If[A\[Equal]2,fxpy[r,\[Theta]],If[A\[Equal]3,fxmy[r,\[Theta]],If[A\[Equal]4,fxyp[r,\[Theta]],fxym[r,\[Theta]]]]]];returnvalue]*)


(*Calculating Source Matrix *)

(*Here, I'm calculating the Source Matrix on r-\[Theta] plane and then making it to one array. 
This source matrix includes the effective source and also it takes into consideration the "jump" between full filed mode \[CapitalPsi]^m and residual field mode Subsuperscript[\[CapitalPsi], R, m] across the worldtube. *)

SourceMatrix[n_, m_, a_, r0_, rminapp_, rmaxapp_, d_]:=Module[{rplus,rminus,rtorstar,r0star,rsmin,rsmax,guess,rstartor,listr,listrsource,\[Theta]smin,\[Theta]smax,imax,jmax,rsourcestart,rsourceend,listSeff1,SourceMatrix1,i,j,r,LeftBC,RightBC,TopBC,BottomBC,SecondTopBC,SecondBottomBC,SecondLeftBC,SourceBC1,SecondRightBC,Sourcelist,\[CapitalDelta]\[Theta],\[Theta]sourcegrid,l,\[CapitalDelta]rstar,rmatrixfactorminus,rmatrixfactorplus,rmin,rmax,\[Theta]start,\[Theta]end},\[CapitalDelta]\[Theta]=\[Pi]/(5(2n+1));\[Theta]sourcegrid=2n+1;l=3n+1;
\[CapitalDelta]rstar=d/(2l+1);
rplus=M+Sqrt[M^2-a^2];rminus=M-Sqrt[M^2-a^2];rtorstar[r_]=r+(2M )/(rplus - rminus) (rplus Log[(r-rplus)/(2M)]-rminus Log[(r-rminus)/(2M)]);
r0star=rtorstar[r0];
rsmin=r0star-\[CapitalDelta]rstar/2-l*\[CapitalDelta]rstar;rsmax=r0star+\[CapitalDelta]rstar/2+l*\[CapitalDelta]rstar;
\[Theta]smin=\[Pi]/2-\[CapitalDelta]\[Theta]/2-n*\[CapitalDelta]\[Theta];\[Theta]smax=\[Pi]/2+\[CapitalDelta]\[Theta]/2+n*\[CapitalDelta]\[Theta];
rmatrixfactorminus=Floor[(rsmin-rminapp)/d]+1;
rmatrixfactorplus=Floor[(rmaxapp-rsmax)/d]+1;
rmin=rsmin-(rmatrixfactorminus-1)*d;
rmax=(rmatrixfactorplus-1)*d+rsmax;
\[Theta]start=\[Theta]smin/\[CapitalDelta]\[Theta];\[Theta]end=\[Theta]smax/\[CapitalDelta]\[Theta];

imax=Round[(rmax-rmin)/\[CapitalDelta]rstar+1];jmax=Floor[\[Pi]/\[CapitalDelta]\[Theta]+1];
guess[rstar_]:=If[rstar<= -2,rplus+2*(1-a^2)^(1/rplus-1/2) Exp[(a^2-rplus+(rplus-1)*rstar)/rplus],If[-2<rstar<= 1000,rplus+2(ProductLog[Exp[1/2 (rstar-rplus)]]),rstar+rplus]];
rstartor[rstar_]:=If[rstar<-50,guess[rstar],r/.FindRoot[rtorstar[r]-rstar==0,{r,guess[rstar]},Method->"Newton"]];
listr=Table[Re[rstartor[rstar]],{rstar,rmin,rmax,\[CapitalDelta]rstar}]//Quiet;listrsource=Table[Re[rstartor[rstar]],{rstar,rsmin,rsmax,\[CapitalDelta]rstar}]//Quiet;
rsourcestart=Floor[(rsmin-rmin)/\[CapitalDelta]rstar+1];rsourceend=rsourcestart+Length[listrsource]-1;
listSeff1=Table[scaledSeffm[m,listrsource[[i]],N[j]],{i,1,Length[listrsource]},{j,\[Theta]smin,\[Theta]smax,\[CapitalDelta]\[Theta]}];

SourceMatrix1=SparseArray[Flatten[Table[{i+rsourcestart-1,j+\[Theta]start}->listSeff1[[i,j]],{i,1,Length[listrsource]},{j,1,\[Theta]end-\[Theta]start+1}]],{imax,jmax}];LeftBC:=SparseArray[Flatten[Table[{i2+rsourcestart-1,\[Theta]start}-> -getfxyCoeff[listrsource[[i2]],N[\[Theta]smin-\[CapitalDelta]\[Theta]],4,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m] listrsource[[i2]] scaled\[CapitalPhi]Pm[m,listrsource[[i2]],N[\[Theta]smin]],{i2,1,Length[listrsource]}]],{imax,jmax}];SecondLeftBC:=SparseArray[Flatten[Table[{i2+rsourcestart-1,\[Theta]start+1}-> getfxyCoeff[listrsource[[i2]],N[\[Theta]smin],5,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m] listrsource[[i2]] scaled\[CapitalPhi]Pm[m,listrsource[[i2]],N[\[Theta]smin-\[CapitalDelta]\[Theta]]],{i2,1,Length[listrsource]}]],{imax,jmax}];RightBC:=SparseArray[Flatten[Table[{i1+rsourcestart-1,\[Theta]end+2}-> -getfxyCoeff [listrsource[[i1]],N[\[Theta]smax+\[CapitalDelta]\[Theta]],5,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m]listrsource[[i1]] scaled\[CapitalPhi]Pm[m,listrsource[[i1]],N[\[Theta]smax]],{i1,1,Length[listrsource]}]],{imax,jmax}];
SecondRightBC:=SparseArray[Flatten[Table[{i1+rsourcestart-1,\[Theta]end+1}-> getfxyCoeff[listrsource[[i1]],N[\[Theta]smax],4,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m] listrsource[[i1]] scaled\[CapitalPhi]Pm[m,listrsource[[i1]],N[\[Theta]smax+\[CapitalDelta]\[Theta]]],{i1,1,Length[listrsource]}]],{imax,jmax}];TopBC:=SparseArray[Flatten[Table[{rsourcestart-1,\[Theta]start+j}-> -getfxyCoeff [listr[[rsourcestart-1]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)],2,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m]listr[[rsourcestart]] scaled\[CapitalPhi]Pm[m,listr[[rsourcestart]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)]],{j,1,\[Theta]end-\[Theta]start+1}]],{imax,jmax}];
SecondTopBC:=SparseArray[Flatten[Table[{rsourcestart,\[Theta]start+j}-> getfxyCoeff[listr[[rsourcestart]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)],3,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m] listr[[rsourcestart-1]] scaled\[CapitalPhi]Pm[m,listr[[rsourcestart-1]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)]],{j,1,\[Theta]end-\[Theta]start+1}]],{imax,jmax}];BottomBC:=SparseArray[Flatten[Table[{rsourceend+1,\[Theta]start+j}-> -getfxyCoeff[listr[[rsourceend+1]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)],3,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m] listr[[rsourceend]] scaled\[CapitalPhi]Pm[m,listr[[rsourceend]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)]],{j,1,\[Theta]end-\[Theta]start+1}]],{imax,jmax}];SecondBottomBC:=SparseArray[Flatten[Table[{rsourceend,\[Theta]start+j}-> getfxyCoeff[listr[[rsourceend]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)],2,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m] listr[[rsourceend+1]] scaled\[CapitalPhi]Pm[m,listr[[rsourceend+1]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)]],{j,1,\[Theta]end-\[Theta]start+1}]],{imax,jmax}];
SourceBC1=LeftBC+SourceMatrix1+SecondLeftBC+RightBC+SecondRightBC+TopBC+SecondTopBC+BottomBC+SecondBottomBC;
Sourcelist=Flatten[ArrayReshape[SourceBC1,{1,imax*jmax}]];
Sourcelist]



(*Calculating  Coupling Matrix *)

(*Here, I'm calculating the Coupling Matrix 
Boundary condition at Subsuperscript[r, max, *] has the error of (1/Subsuperscript[r, max, *])^4  e^(i\[Omega] Subsuperscript[r, max, *]) *)

CouplingMatrix[n_, m_, a_, r0_, rminapp_, rmaxapp_, d_] :=Module[{rplus,rminus,rtorstar,r0star,rstarmin1,rstarmax1,rsmin1,rsmax1,\[Theta]smin,\[Theta]smax,imax,jmax,listr1,listrsource1,MxBC,rstartor,listr,listrsource,guess,\[Omega],rsmin,rsmax,rmin,rmax,\[CapitalDelta]\[Theta],\[Theta]sourcegrid,l,\[CapitalDelta]rstar,rsourcegrid,rmatrixfactorminus,rmatrixfactorplus,listTheta, diag, rightDiag, leftDiag,
     rightSkipDiag, leftSkipDiag,xBC1,xBC2,xBC3,yBC1,yBC2,yBC3,xBC1min,xBC1max,xBC2min,xBC2max,xBC3min,xBC3max,xBC4max,xBC5max,listr2},\[Omega]=m*\[CapitalOmega];\[CapitalDelta]\[Theta]=\[Pi]/(5(2n+1));\[Theta]sourcegrid=2n+1;l=3n+1;\[CapitalDelta]rstar=d/(2l+1);
rsourcegrid=2l+1;
rplus=M+Sqrt[M^2-a^2];rminus=M-Sqrt[M^2-a^2];rtorstar[r_]=r+(2M )/(rplus - rminus) (rplus Log[(r-rplus)/(2M)]-rminus Log[(r-rminus)/(2M)]);
r0star=rtorstar[r0];
rsmin=r0star-\[CapitalDelta]rstar/2-l*\[CapitalDelta]rstar;rsmax=r0star+\[CapitalDelta]rstar/2+l*\[CapitalDelta]rstar;
\[Theta]smin=\[Pi]/2-\[CapitalDelta]\[Theta]/2-n*\[CapitalDelta]\[Theta];\[Theta]smax=\[Pi]/2+\[CapitalDelta]\[Theta]/2+n*\[CapitalDelta]\[Theta];
rmatrixfactorminus=Floor[(rsmin-rminapp)/d]+1;
rmatrixfactorplus=Floor[(rmaxapp-rsmax)/d]+1;
rmin=rsmin-(rmatrixfactorminus-1)*d;
rmax=(rmatrixfactorplus-1)*d+rsmax;

imax=Round[(rmax-rmin)/\[CapitalDelta]rstar];jmax=Floor[\[Pi]/\[CapitalDelta]\[Theta]];
guess[rstar_]:=If[rstar<= -2,rplus+2*(1-a^2)^(1/rplus-1/2) Exp[(a^2-rplus+(rplus-1)*rstar)/rplus],If[-2<rstar<= 1000,rplus+2(ProductLog[Exp[1/2 (rstar-rplus)]]),rstar+rplus]];
rstartor[rstar_]:=If[rstar<-50,guess[rstar],r/.FindRoot[rtorstar[r]-rstar==0,{r,guess[rstar]},Method->"Newton"]];
listr1=Table[Re[rstartor[rstar]],{rstar,rmin,rmax,\[CapitalDelta]rstar}][[2;;-2]]//Quiet;


listrsource=Table[Re[rstartor[rstar]],{rstar,rsmin,rsmax,\[CapitalDelta]rstar}]//Quiet;
 listTheta = Table[j * \[CapitalDelta]\[Theta], {j, 1, jmax - 1}];
   
diag = Flatten[ArrayPad[Transpose[Table[(a^2 m^2 listr1^2 (a^2 + 
      listr1 (-2 M + listr1)) \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 \[CapitalOmega]^2 + (2 listr1^2 (a^2 + listr1
       (-2 M + listr1)) \[CapitalDelta]rstar^2 + 2 (listr1^2 (a^2 + listr1^2)^2 - (a^2 - 
      I a m listr1 - M listr1) (a^2 + listr1 (-2 M + listr1)) \[CapitalDelta]rstar^2) \[CapitalDelta]\[Theta]^
      2 + 4 a m^2 M listr1^3 \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 \[CapitalOmega] - m^2 listr1^2 (a^2 + listr1^2
      )^2 \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 \[CapitalOmega]^2) Csc[theta]^2 + m^2 listr1^2 (a^2 + listr1 (-2 
      M + listr1)) \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 Csc[theta]^4) / (listr1^2 \[CapitalDelta]rstar^2 \[CapitalDelta]\[Theta]^2 (-
      a^2 (a^2 + listr1 (-2 M + listr1)) + (a^2 + listr1^2)^2 Csc[theta]^2)
      ), {theta, listTheta}]], 1]];
    rightDiag = Flatten[ArrayPad[Transpose[Table[(a^2 + listr1 (-2 M 
      + listr1)) (2 + \[CapitalDelta]\[Theta] Cot[theta]) / (2 \[CapitalDelta]\[Theta]^2 (-(a^2 + listr1^2)^2 + a^2 (
      a^2 + listr1 (-2 M + listr1)) Sin[theta]^2)), {theta, listTheta}]], 1
      ]];
    leftDiag = Flatten[ArrayPad[Transpose[Table[-(a^2 + listr1 (-2 M 
      + listr1)) (-2 + \[CapitalDelta]\[Theta] Cot[theta]) / (2 \[CapitalDelta]\[Theta]^2 (-(a^2 + listr1^2)^2 + a^2 
      (a^2 + listr1 (-2 M + listr1)) Sin[theta]^2)), {theta, listTheta}]], 
      1]];
    rightSkipDiag = Flatten[ArrayPad[Transpose[Table[(-listr1 (a^2 + 
      listr1^2)^2 + a (a^3 - I a^2 m listr1 - I m listr1^3 + a listr1 (-2 M
       + listr1)) \[CapitalDelta]rstar) / (listr1 \[CapitalDelta]rstar^2 ((a^2 + listr1^2)^2 - a^2 (a^2
       + listr1 (-2 M + listr1)) Sin[theta]^2)), {theta, listTheta}]], {{1,
       0}, {1, 1}}]];
    leftSkipDiag = Flatten[ArrayPad[Transpose[Table[-(listr1 (a^2 + listr1
      ^2)^2 + a (a^3 - I a^2 m listr1 - I m listr1^3 + a listr1 (-2 M + listr1
      )) \[CapitalDelta]rstar) / (listr1 \[CapitalDelta]rstar^2 ((a^2 + listr1^2)^2 - a^2 (a^2 + listr1
       (-2 M + listr1)) Sin[theta]^2)), {theta, listTheta}]], {{0, 1}, {1, 
      1}}]];
xBC1min=3-2I \[CapitalDelta]rstar \[Omega];xBC2min=-4;xBC3min=1;

xBC1max=5-12 I \[CapitalDelta]rstar \[Omega]-9 \[CapitalDelta]rstar^2 \[Omega]^2+2 I \[CapitalDelta]rstar^3 \[Omega]^3;
xBC2max=-18 + 30 I \[CapitalDelta]rstar \[Omega] + 12 \[CapitalDelta]rstar^2 \[Omega]^2;
xBC3max=24-24 I \[CapitalDelta]rstar \[Omega] - 3 \[CapitalDelta]rstar^2 \[Omega]^2;
xBC4max=-14+6 I \[CapitalDelta]rstar \[Omega];
xBC5max=3;
yBC1=If[m==0,3,1];yBC2=If[m==0,-4,0];yBC3=If[m==0,1,0];
MxBC=SparseArray[{Band[{1, 1}, {(jmax + 1), (jmax + 1)}] -> xBC1min,
Band[{1,  jmax + 2}, {(jmax + 1), 2 (1+jmax)}] -> xBC2min,
Band[{1,3+2 jmax}, {(jmax + 1), 3 (1+jmax)}] -> xBC3min,
Band[{1+imax+imax jmax,  (imax + 1) (jmax + 1)-jmax}, {(imax + 1) (jmax + 1), (imax + 1) (jmax + 1)}] -> xBC1max,
Band[{1+imax+imax jmax, (imax + 1) (jmax + 1)-jmax-1-jmax}, {(imax + 1) (jmax + 1), (imax + 1) (jmax + 1)-jmax-1}] -> xBC2max,
Band[{1+imax+imax jmax,(imax + 1) (jmax + 1)-jmax-1-jmax-1-jmax}, {(imax + 1) (jmax + 1),(imax + 1) (jmax + 1)-jmax-1-jmax-1}] -> xBC3max,
Band[{1+imax+imax jmax,(imax + 1) (jmax + 1)-jmax-1-jmax-1-jmax-1-jmax},{(imax + 1) (jmax + 1),(imax + 1) (jmax + 1)-jmax-1-jmax-1-jmax-1}] ->xBC4max,
Band[{1+imax+imax jmax,(imax + 1) (jmax + 1)-jmax-1-jmax-1-jmax-1-jmax-1-jmax},{(imax + 1) (jmax + 1),(imax + 1) (jmax + 1)-jmax-1-jmax-1-jmax-1-jmax-1}] ->xBC5max,
Band[{2+jmax,2+jmax},{imax+(-1+imax) jmax,imax+(-1+imax) jmax},1+jmax]->yBC1,Band[{2(1+jmax),2(1+jmax)},{(1+imax-1)(1+jmax),(1+imax-1)(1+jmax)},1+jmax]->yBC1,Band[{2+jmax,3+jmax},{imax+(-1+imax) jmax,imax+(-1+imax) jmax+1},1+jmax]->yBC2,Band[{2(1+jmax),2(1+jmax)-1},{(1+imax-1)(1+jmax),(1+imax-1)(1+jmax)-1},1+jmax]->yBC2,Band[{2+jmax,4+jmax},{imax+(-1+imax) jmax,imax+(-1+imax) jmax+2},1+jmax]->yBC3,Band[{2(1+jmax),2(1+jmax)-2},{(1+imax-1)(1+jmax),(1+imax-1)(1+jmax)-2},1+jmax]->yBC3}, {(imax + 1) (jmax + 1), (imax + 1) (jmax + 1)}, 0] + DiagonalMatrix[SparseArray[
      {i_} :> diag[[i]], {Length[diag]}, 0], 0] + DiagonalMatrix[SparseArray[
      {i_} :> rightDiag[[i]], {Length[rightDiag] - 1}, 0], 1]+ DiagonalMatrix[
      SparseArray[{i_} :> leftDiag[[i + 1]], {Length[leftDiag] - 1}, 0], -1
      ] + DiagonalMatrix[SparseArray[{i_} :> rightSkipDiag[[i]], {Length[rightSkipDiag
      ]}, 0], jmax + 1] + DiagonalMatrix[SparseArray[{i_} :> leftSkipDiag[[
      i]], {Length[leftSkipDiag]}, 0], -jmax - 1];
MxBC]
      


(*Linear Solve *)

(*Here, I'm solving the linear system, Subscript[M, coupling] * Subscript[M, \[CapitalPsi]^m] = Subscript[M, source^m]*)

\[Psi]FDA[n_, m_, a_, r0_, rminapp_, rmaxapp_, d_]:=Module[{final\[Psi]1,jmax,\[CapitalDelta]\[Theta]},\[CapitalDelta]\[Theta]=\[Pi]/(5(2n+1));jmax=Floor[\[Pi]/\[CapitalDelta]\[Theta]+1];

final\[Psi]1:=Partition[LinearSolve[CouplingMatrix[n,m,a,r0,rminapp,rmaxapp,d],SourceMatrix[n,m,a,r0,rminapp,rmaxapp,d],Method->"Pardiso"],jmax];
final\[Psi]1];
\[Psi]sol=\[Psi]FDA[n,m,a,r0,rminapp,rmaxapp,d]


(*Reassigning i-j grids to r-\[Theta] grids *)

(*Here, I'm re-assigning \[Psi][i, j] to \[Psi][(r^*),\[Theta]] 
The solution includes Subscript[i, max], Subscript[j, max], \[CapitalDelta]r^*, \[CapitalDelta]\[Theta], which will be used later . *)

\[Psi]atr\[Theta][n_, m_, a_, r0_, rminapp_, rmaxapp_, d_]:=Module[{imax,rstarmax,rstarmin,r0star,jmax,rplus,rminus,rtorstar,listr,rstartor,listrstar,list\[Theta],guess,function\[CapitalPsi],d\[CapitalPsi]r,\[CurlyPhi],\[CapitalDelta],\[Psi]inr\[Theta]1try1,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],rmatrixfactorminus,rmatrixfactorplus,rmin,rsmin,rsmax,rmax,l},\[CapitalDelta]\[Theta]=\[Pi]/(5(2n+1));l=3n+1;\[CapitalDelta]rstar=d/(2l+1);
rplus=M+Sqrt[M^2-a^2];rminus=M-Sqrt[M^2-a^2];rtorstar[r_]=r+(2M )/(rplus - rminus) (rplus Log[(r-rplus)/(2M)]-rminus Log[(r-rminus)/(2M)]);
r0star=rtorstar[r0];rsmin=r0star-\[CapitalDelta]rstar/2-l*\[CapitalDelta]rstar;rsmax=r0star+\[CapitalDelta]rstar/2+l*\[CapitalDelta]rstar;
rmatrixfactorminus=Floor[(rsmin-rminapp)/d]+1;
rmatrixfactorplus=Floor[(rmaxapp-rsmax)/d]+1;
rmin=rsmin-(rmatrixfactorminus-1)*d;
rmax=(rmatrixfactorplus-1)*d+rsmax;
guess[rstar_]:=If[rstar<= -2,rplus+2*(1-a^2)^(1/rplus-1/2) Exp[(a^2-rplus+(rplus-1)*rstar)/rplus],If[-2<rstar<= 1000,rplus+2(ProductLog[Exp[1/2 (rstar-rplus)]]),rstar+rplus]];
rstartor[rstar_]:=If[rstar<-50,guess[rstar],r/.FindRoot[rtorstar[r]-rstar==0,{r,guess[rstar]},Method->"Newton"]];
imax=Round[(rmax-rmin)/\[CapitalDelta]rstar+1];jmax=Floor[\[Pi]/\[CapitalDelta]\[Theta]+1];listr=Table[Re[rstartor[rstar]],{rstar,rmin,rmax,\[CapitalDelta]rstar}]//Quiet;
listrstar=Table[rtorstar[listr[[i]]],{i,1,Length[listr]}];
list\[Theta]=Table[(j-1)\[CapitalDelta]\[Theta],{j,1,jmax}];
\[CurlyPhi][r_]=a/(rplus-rminus) Log[(r-rplus)/(r-rminus)]; 
\[CapitalDelta][r_]=r^2-2 M r+a^2;
{\[Psi]inr\[Theta]1try1=Flatten[Table[{listrstar[[i]],list\[Theta][[j]],2*\[Psi]sol[[i,j]]},{i,1,imax},{j,1,Length[list\[Theta]]}],1],imax,jmax,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta]}]


(*Calculating Scalar Self Force*)

SSF[n_, m_, a_, r0_, rminapp_, rmaxapp_, d_]:=Module[{rplus,rminus,\[CurlyPhi],\[Phi],\[CapitalDelta],\[Psi]atr\[Theta]sol,rtorstar,r0star,imax,jmax,\[CapitalDelta]rstar1,\[CapitalDelta]\[Theta]1,function\[Psi],d\[CapitalPsi]r,d\[Phi]r1,d\[Phi]\[Theta],d\[Phi]\[Phi]},rplus=M+Sqrt[M^2-a^2];rminus=M-Sqrt[M^2-a^2];
\[CurlyPhi][r_]:=\[Phi]+a/(rplus-rminus) Log[(r-rplus)/(r-rminus)]; 
\[Phi]=0; t=0;
\[CapitalDelta][r_]:=r^2-2 M r+a^2;
\[Psi]atr\[Theta]sol=\[Psi]atr\[Theta][n,m,a,r0,rminapp,rmaxapp,d];
rtorstar[r_]=r+(2M )/(rplus - rminus) (rplus Log[(r-rplus)/(2M)]-rminus Log[(r-rminus)/(2M)]);
r0star=rtorstar[r0]; 
imax=\[Psi]atr\[Theta]sol[[2]];
jmax=\[Psi]atr\[Theta]sol[[3]];(*jmax is always even number*)
\[CapitalDelta]rstar1=\[Psi]atr\[Theta]sol[[4]];
\[CapitalDelta]\[Theta]1=\[Psi]atr\[Theta]sol[[5]];
function\[Psi]=Interpolation[\[Psi]atr\[Theta][n,m,a,r0,rminapp,rmaxapp,d][[1]]];
d\[CapitalPsi]r=(-function\[Psi][r0star+2\[CapitalDelta]rstar1,\[Pi]/2]+8function\[Psi][r0star+\[CapitalDelta]rstar1,\[Pi]/2]-8function\[Psi][r0star-\[CapitalDelta]rstar1,\[Pi]/2]+function\[Psi][r0star-2\[CapitalDelta]rstar1,\[Pi]/2])/(12\[CapitalDelta]rstar1);
d\[Phi]r1=(1/r0 function\[Psi][r0star,\[Pi]/2] I m a/((r0-rminus)(r0-rplus))+1/r0 *(r0^2+a^2)/\[CapitalDelta][r0] d\[CapitalPsi]r - function\[Psi][r0star,\[Pi]/2]  *1/r0^2);
d\[Phi]\[Theta]=(-function\[Psi][r0star,\[Pi]/2+2\[CapitalDelta]\[Theta]1]+8function\[Psi][r0star,\[Pi]/2+\[CapitalDelta]\[Theta]1]-8function\[Psi][r0star,\[Pi]/2-\[CapitalDelta]\[Theta]1]+function\[Psi][r0star,\[Pi]/2-2\[CapitalDelta]\[Theta]1])/(12\[CapitalDelta]\[Theta]1);
d\[Phi]\[Phi]=1/r0 function\[Psi][r0star,\[Pi]/2]*(I m);
Return[{d\[Phi]r1,d\[Phi]\[Phi]}]]



(*You have to be careful when a \[NotEqual]0 
d\[Phi]r1= (1/r0 function\[Psi][r0star,\[Pi]/2] I m a/((r0-rminus)(r0-rplus))+1/r0 *(r0^2+a^2)/\[CapitalDelta][r0] d\[CapitalPsi]r - function\[Psi][r0star,\[Pi]/2]  *1/r0^2)  * e^(i m \[CurlyPhi])
\[CurlyPhi] = \[Phi] + \[CapitalDelta]\[Phi][r] 
\[CapitalDelta]\[Phi][r] = a/(r^+-r^-) ln|(r-r^+)/(r-r^-)| 
and therefore, when a \[NotEqual] 0, \[CurlyPhi] \[NotEqual] 0, and therefore e^im\[CurlyPhi] \[NotEqual] 1 *)


(*Scalar Self Force Results *)

ssf=SSF[n,m,a,r0,rminapp,rmaxapp,d] ; 
Fr = ssf[[1]]
F\[Phi] = ssf[[2]]
