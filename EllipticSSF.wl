(* ::Package:: *)

(*Calculating Scalar Self Force Via Elliptic PDEs *)
(*This is a  Mathematica Notebook to calculate the scalar self force via Elliptic PDEs. *)


(*You can change those values to try different mode m, spin parameter a, mass M, and resolutions. *)

m=; (* mode m *)
n=; (* resolution *)
rminapp=;(*approximate value for the r^* minimum *)
rmaxapp=;(*approximate value for the r^* maximum  *)
a=; (* spin parameter *)
r0=;(* location of the worldline *)
q=1; (* sclar charge *)
d=3.6; (*width of the worldtube *)
M=1; (* mass of the central black hole *)
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

SourceCombinedNew[n_, m_, a_, r0_, rminapp_, rmaxapp_, d_]:=Module[{rplus,rminus,rtorstar,r0star,rsmin,rsmax,guess,rstartor,listr,listrsource,\[Theta]smin,\[Theta]smax,imax,jmax,rsourcestart,rsourceend,listSeff1,SourceMatrix1,i,j,r,LeftBC,RightBC,TopBC,BottomBC,SecondTopBC,SecondBottomBC,SecondLeftBC,SourceBC1,SecondRightBC,Sourcelist,\[CapitalDelta]\[Theta],\[Theta]sourcegrid,l,\[CapitalDelta]rstar,rmatrixfactorminus,rmatrixfactorplus,rmin,rmax,\[Theta]start,\[Theta]end},\[CapitalDelta]\[Theta]=\[Pi]/(5(2n+1));\[Theta]sourcegrid=2n+1;l=3n+1;
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

SourceMatrix1=SparseArray[Flatten[Table[{i+rsourcestart-1,j+\[Theta]start}->listSeff1[[i,j]],{i,1,Length[listrsource]},{j,1,\[Theta]end-\[Theta]start+1}]],{imax,jmax}];LeftBC:=SparseArray[Flatten[Table[{i2+rsourcestart-1,\[Theta]start}-> -getfxyNew[listrsource[[i2]],N[\[Theta]smin-\[CapitalDelta]\[Theta]],4,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m] listrsource[[i2]] scaled\[CapitalPhi]Pm[m,listrsource[[i2]],N[\[Theta]smin]],{i2,1,Length[listrsource]}]],{imax,jmax}];SecondLeftBC:=SparseArray[Flatten[Table[{i2+rsourcestart-1,\[Theta]start+1}-> getfxyNew[listrsource[[i2]],N[\[Theta]smin],5,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m] listrsource[[i2]] scaled\[CapitalPhi]Pm[m,listrsource[[i2]],N[\[Theta]smin-\[CapitalDelta]\[Theta]]],{i2,1,Length[listrsource]}]],{imax,jmax}];RightBC:=SparseArray[Flatten[Table[{i1+rsourcestart-1,\[Theta]end+2}-> -getfxyNew [listrsource[[i1]],N[\[Theta]smax+\[CapitalDelta]\[Theta]],5,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m]listrsource[[i1]] scaled\[CapitalPhi]Pm[m,listrsource[[i1]],N[\[Theta]smax]],{i1,1,Length[listrsource]}]],{imax,jmax}];
SecondRightBC:=SparseArray[Flatten[Table[{i1+rsourcestart-1,\[Theta]end+1}-> getfxyNew[listrsource[[i1]],N[\[Theta]smax],4,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m] listrsource[[i1]] scaled\[CapitalPhi]Pm[m,listrsource[[i1]],N[\[Theta]smax+\[CapitalDelta]\[Theta]]],{i1,1,Length[listrsource]}]],{imax,jmax}];TopBC:=SparseArray[Flatten[Table[{rsourcestart-1,\[Theta]start+j}-> -getfxyNew [listr[[rsourcestart-1]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)],2,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m]listr[[rsourcestart]] scaled\[CapitalPhi]Pm[m,listr[[rsourcestart]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)]],{j,1,\[Theta]end-\[Theta]start+1}]],{imax,jmax}];
SecondTopBC:=SparseArray[Flatten[Table[{rsourcestart,\[Theta]start+j}-> getfxyNew[listr[[rsourcestart]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)],3,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m] listr[[rsourcestart-1]] scaled\[CapitalPhi]Pm[m,listr[[rsourcestart-1]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)]],{j,1,\[Theta]end-\[Theta]start+1}]],{imax,jmax}];BottomBC:=SparseArray[Flatten[Table[{rsourceend+1,\[Theta]start+j}-> -getfxyNew[listr[[rsourceend+1]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)],3,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m] listr[[rsourceend]] scaled\[CapitalPhi]Pm[m,listr[[rsourceend]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)]],{j,1,\[Theta]end-\[Theta]start+1}]],{imax,jmax}];SecondBottomBC:=SparseArray[Flatten[Table[{rsourceend,\[Theta]start+j}-> getfxyNew[listr[[rsourceend]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)],2,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],m] listr[[rsourceend+1]] scaled\[CapitalPhi]Pm[m,listr[[rsourceend+1]],N[\[Theta]smin+\[CapitalDelta]\[Theta]*(j-1)]],{j,1,\[Theta]end-\[Theta]start+1}]],{imax,jmax}];
SourceBC1=LeftBC+SourceMatrix1+SecondLeftBC+RightBC+SecondRightBC+TopBC+SecondTopBC+BottomBC+SecondBottomBC;
Sourcelist=Flatten[ArrayReshape[SourceBC1,{1,imax*jmax}]];
Sourcelist]

SM=SourceCombinedNew[n, m, a, r0, rminapp, rmaxapp, d;]


GetParam[n_, m_, a_, r0_, rminapp_, rmaxapp_, d_] :=Module[{rplus,rminus,rtorstar,r0star,rstarmin1,rstarmax1,rsmin1,rsmax1,\[Theta]smin,\[Theta]smax,imax,jmax,listr1,listrsource1,MxBC,rstartor,listr,listrsource,guess,\[Omega],rsmin,rsmax,rmin,rmax,\[CapitalDelta]\[Theta],\[Theta]sourcegrid,l,\[CapitalDelta]rstar,rsourcegrid,rmatrixfactorminus,rmatrixfactorplus,listTheta,list\[Theta],listrstar},\[Omega]=m*\[CapitalOmega];\[CapitalDelta]\[Theta]=\[Pi]/(5(2n+1));\[Theta]sourcegrid=2n+1;l=3n+1;\[CapitalDelta]rstar=d/(2l+1);
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

listr=Table[Re[rstartor[rstar]],{rstar,rmin,rmax,\[CapitalDelta]rstar}]//Quiet;
listrstar=Table[rtorstar[listr[[i]]],{i,1,Length[listr]}];
listr1=listr[[2;;-2]]//Quiet;
listrsource=Table[Re[rstartor[rstar]],{rstar,rsmin,rsmax,\[CapitalDelta]rstar}]//Quiet;list\[Theta]=Table[(j-1)\[CapitalDelta]\[Theta],{j,1,jmax+1}];
listTheta=list\[Theta][[2;;-2]];{\[Omega],\[CapitalDelta]\[Theta],\[CapitalDelta]rstar,imax,jmax,listr1,listrsource,listTheta,list\[Theta],listrstar,r0star}]


CouplingMatrixNew[n_, m_, a_, r0_, rminapp_, rmaxapp_, d_] :=Module[{\[Omega],\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],listTheta,imax,jmax,listr1,diag,rightDiag,leftDiag,rightSkipDiag,leftSkipDiag,xBC1min,xBC2min,xBC3min,xBC1max,xBC2max,xBC3max,xBC4max,xBC5max,yBC1,yBC2,yBC3,MxBC},
\[Omega]=GetParam[n,m,a,r0,rminapp,rmaxapp,d][[1]];
\[CapitalDelta]rstar=GetParam[n,m,a,r0,rminapp,rmaxapp,d][[3]];
\[CapitalDelta]\[Theta]=GetParam[n,m,a,r0,rminapp,rmaxapp,d][[2]];
listTheta=GetParam[n,m,a,r0,rminapp,rmaxapp,d][[8]];
imax=GetParam[n,m,a,r0,rminapp,rmaxapp,d][[4]];
jmax=GetParam[n,m,a,r0,rminapp,rmaxapp,d][[5]];
listr1=GetParam[n,m,a,r0,rminapp,rmaxapp,d][[6]];
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
      


SSF[n_, m_, a_, r0_, rminapp_, rmaxapp_, d_]:=Module[{rplus,rminus,\[CurlyPhi],\[Phi],\[CapitalDelta],\[Psi]valuesatr\[Theta]New,rtorstar,r0star,imax,jmax,\[CapitalDelta]rstar,\[CapitalDelta]\[Theta],function\[Psi],d\[CapitalPsi]r,d\[Phi]r1,d\[Phi]\[Theta],d\[Phi]\[Phi]},
rplus=M+Sqrt[M^2-a^2];rminus=M-Sqrt[M^2-a^2];
\[CapitalDelta]\[Theta]=GetParam[n,m,a,r0,rminapp,rmaxapp,d][[2]];
\[CapitalDelta]rstar=GetParam[n,m,a,r0,rminapp,rmaxapp,d][[3]];
r0star=GetParam[n,m,a,r0,rminapp,rmaxapp,d][[11]];
\[CapitalDelta][r_]:=r^2-2 M r+a^2;
\[CurlyPhi][r_]:=\[Phi]+a/(rplus-rminus) Log[(r-rplus)/(r-rminus)]; 
\[Phi]=0; 

function\[Psi]=Interpolation[\[Psi]atr\[Theta]New[n,m,a,r0,rminapp,rmaxapp,d][[1]]];
d\[CapitalPsi]r=(-function\[Psi][r0star+2\[CapitalDelta]rstar,\[Pi]/2]+8function\[Psi][r0star+\[CapitalDelta]rstar,\[Pi]/2]-8function\[Psi][r0star-\[CapitalDelta]rstar,\[Pi]/2]+function\[Psi][r0star-2\[CapitalDelta]rstar,\[Pi]/2])/(12\[CapitalDelta]rstar);
d\[Phi]r1=(1/r0 function\[Psi][r0star,\[Pi]/2] I m a/((r0-rminus)(r0-rplus))+1/r0 *(r0^2+a^2)/\[CapitalDelta][r0] d\[CapitalPsi]r - function\[Psi][r0star,\[Pi]/2]  *1/r0^2)*Exp[I m \[CurlyPhi][r0]];
d\[Phi]\[Theta]=(-function\[Psi][r0star,\[Pi]/2+2\[CapitalDelta]\[Theta]]+8function\[Psi][r0star,\[Pi]/2+\[CapitalDelta]\[Theta]]-8function\[Psi][r0star,\[Pi]/2-\[CapitalDelta]\[Theta]]+function\[Psi][r0star,\[Pi]/2-2\[CapitalDelta]\[Theta]])/(12\[CapitalDelta]\[Theta]);
d\[Phi]\[Phi]=1/r0 function\[Psi][r0star,\[Pi]/2]*(I m)*Exp[I m \[CurlyPhi][r0]];
Return[{d\[Phi]r1,d\[Phi]\[Phi]}]]



ssf=SSF[n,m,a,r0,rminapp,rmaxapp,d] ; 
SSFvalue=If[m==0,ssf,2*ssf]
