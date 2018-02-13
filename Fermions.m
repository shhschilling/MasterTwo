

(* Gamma matrices are enoted by Dirac[p] or Dirac[mu], 
   the unit matrix by Unit, Gamma5 by Gamma5  *)

(* Define the right R and the left projector L to be taken to the left *)

(* The Scalar product is denoted by Scal[p,q]; its output is p\[CenterDot]q \
*)

(* This version works with replacement rules only, hoping for more
   speed .. *)

(* Momenta, Indices and Masses have to be declared with
   DeclareMomentum[p],... *)



BeginPackage["Fermions`"]





Unprotect["Fermions`*"];



(* Startup message *)

     (*Print["\n ***  \t Fermions \t ***"];*)
     (*Print[" For a general information on this package try 'FermionsInfo[]'"];*)

FermionsInfo[] := (
    Print["\n *** \t Fermions \t ***"];
    Print[" Fermions is a package that can simplify Dirac expressions in d \
dimensions with an anticommuting Gamma5."];
    Print["\n Available commands:"];
    Print[" (For any of these commands '?command' should help.)\n"];
    Print[Names["Fermions`*"]];
		   );








 TraceMoreThanFour::usage = "";


 Dirac::usage =
    "Dirac[p] denotes p-slash; Dirac[mu] is gamma_mu; Products: Dirac[L, p + \
m , mu ] 
     (have to declare p with 'DeclareMomentum[p]', mu with 'DeclareIndex[mu]' \
and 
     m with 'DeclareMass[m]' first).";

 d::usage =
    "d is the space-time dimension. To set it globally, have to unprotect it \
first: 
      Unprotect[d]; d = 4;";

 L::usage =
    "L denotes 1/2(1-Gamma5); usage: Dirac[L]; Projectors are
     anticommutated to the left automatically.";

 R::usage =
    "R denotes 1/2(1+Gamma5); usage: Dirac[R]; Projectors are
     anticommutated to the left automatically.";

 Unit::usage =
    "Unit denotes the unit matrix; usage: Dirac[Unit] or Dirac[]";

 Gamma5::usage =
    "Gamma5 denotes gamma_5; usage: Dirac[Gamma5]. Gamma5 is 
     anticommutated to the left automatically.";

 Scal::usage =
    "Scal[p,q] is the scalar product.";

 ContractIndex::usage =  
    "ContractIndex[expr,{mu,nu,...}] contracts all mu's etc. May have
     to Expand (or use DiracLinearity) first.";  

 DiracSort::usage =
    "DiracSort[expr,reflist] orders any sufficiently simple expression
     of gammas in the order specified in reflist (a list containing
     all the momenta and indices appearing in expr). It may be a good
     idea to use DiracCollect[expr] first. Make sure no Projectors (L,R
     or Gamma5) appear in reflist! Also make sure that all momenta
     appear in reflist.";

 DiracCollect::usage =
    "DiracCollect[expr] collects Dirac[] expressions.
     Can also use DiracCollect[expr,func], then func is applied to the \
coefficients
     e.g. DiracCollect[expr,Simplify].";

 DiracFactor::usage =
    "DiracFactor[expr] collects all Dirac[]'s and factors the coefficient. \
(It is an alias
     for DiracCollect[expr,Factor])";

 DeclareMass::usage =
    "Declare any mass m you wish to use with DeclareMass[m]; Can use it then \
as
     Dirac[m] or Dirac[p + m]. Side effect: Conjugate[m] will evaluate as \
m.";

 DeclareMomentum::usage =
    "Declare any momentum p you wish to use with DeclareMomentum[p]";

 DeclareIndex::usage =
    "Declare any index mu you wish to use with DeclareIndex[mu]. The indices
     are stored in ListOfIndices.";

 DeclarePolarizationVector::usage = 
    "DeclarePolarizationVector[{e,k}] (k has to be declared as a momentum \
first) 
     declares e as a momentum, sets Scal[k,e]=0 and Scal[Conjugate[e],k]=0. 
     DiracAdjunction[e] will return Conjugate[e] (and not e itself as for a \
momentum). 
     DeclarePolarizationVector[e] does the same but sets no scalar products \
to zero.";

 ListOfIndices::usage = 
    "ListOfIndices is the list of all declared Indices.";

 DiracLinearity::usage =
    "DiracLinearity[expr] expands all sums within Dirac[] and takes \
prefactors of masses,
     momenta and indices out of Dirac[]. It does the same for Scal[]";

 DiracAlgebra::usage =
    "DiracAlgebra[expr] does some of the standard Gamma algebra. May have to
     contract a few indices by hand afterwards (if they are too far apart).";

 ContractAllIndices::usage =
    "ContractAllIndices[expr] contracts all silent indices (may have to use 
     'DiracLinearity' first).";

 UseDiracEquation::usage = 
    "UseDiracEquation[{p,mp},expr,{q,mq}]
      sorts expr and uses the Dirac eq (as in \bar{u}(p) expr u(q); not
      v's).  For Antiparticles use
      UseDiracEquation[{p,-mp},expr,{q,-mq}] (as in \bar{v}(p) expr v(q) ).
      In the same spirit: UseDiracEquation[{p,mp},expr,{}] and
      UseDiracEquation[{},expr,{q,mq}] can be used.";

 MakeReflist::usage = "MakeReflist[expr,{left,right}] makes a reference list \
of all
    the momenta, indices and masses in the Dirac expression. left will be on \
the left
    of this list, right on the right (UseDiracEquation uses this function)";

 DiracTrace::usage =
    "DiracTrace[Dirac[...]] represents the trace over the Dirac expression; \
the trace is
     not evaluated. Non-Dirac expressions may taken out of DiracTrace[...] \
using
     DiracTraceLinearity. To evaluate the
     trace, DiracTraceAlgebra has to be used.";

 DiracTraceLinearity::usage =
    "DiracTraceLinearity[DiracTrace[...]] should be obvious...!";

 DiracTraceAlgebra::usage =
    "DiracTraceAlgebra[expr] evaluates the traces over all \
DiracTrace[Dirac[...]] expressions.
     It's a good idea to use DiracAlgebra before tracing an expression.
     To simplifty the expression afterwards use Expand, ContractScal and
     EpsilonSort.
     WARNING: If the endresult contains 'Epsilon', these traces have been
     calculated in d=4 only!";

 DiracTraceList::usage =
    "The list of the basic rules for DiracTraceAlgebra.";

 LearnDiracTraceRule::usage =
    "LearnDiracTraceRule[Dirac[k1,k2,...,kn]] increases the speed of \
calculations of traces with
     projectors or with a larger amount of momenta. 
     This routine adds rules to 'ExtendedDiracTraceList'. Make sure that only \
momenta, indices or 
     projectors are input to this command. If Dirac[k1,...k10] is entered the \
routine will also
     learn the rule for any shorter expression of this form (e.g. \
Dirac[k1,...k8] and Dirac[k1,...,k6]).
     Also, make sure that there are no rules for the scalarproducts of these \
momenta, such as Scal[k1,k2] = 0!!";

 ExtendedDiracTraceList::usage =
    "GamTrace[] knows only few rules. In this list you can supply more rules \
to make GamTrace faster
     (this is highly recommended!!). Check also 'LearnGamTraceRule'";

 Epsilon::usage =
    "Epsilon[a,b,c,d] is the completely antisymmetric tensor in 4 \
dimensions."
    eps::usage="eps is the expansion parameter."
  ContractScalEps::usage =
    "ContractScalEps[expr]  contracts things like Epsilon[a,___] Scal[a,_]."

  EpsilonSort::usage =
    "EpsilonSort[expr,reflist] sorts expressions in Epsilon[___] according
     to reflist."

  DiracAdjunction::usage = 
    "Is what you'd guess!";

 DiracSquare::usage =
    "DiracSquare[expr,one,two] returns the trace of expr * one * \
DiracAdjunction[expr] * two.
     one and two should be Dirac[___] expressions.";

 DiracProduct::usage =
    "DiracProduct[a Dirac[b], c Dirac[d],... ] returns: a c Dirac[b,d,...]";

 MomentumQ::usage =
    "Is what you'd guess...";

 ProjectorQ::usage =
    "R, L, Gamma5 and Unit are treated as projectors. ProjectorQ[R] will \
return 'True'.";

 MassQ::usage =
    "Is what you'd guess...";

 IndexQ::usage =
    "Is what you'd guess...";

 DiracVectorQ::usage =
    "DiracVectorQ = IndexQ || MomentumQ";

 PolarizationVectorQ::usage = 
    "Is what you'd guess...";

 DiracScalExpand::usage =
    "DiracScalExpand[expr] expands all arguments in Dirac[] and Scal[]";

 IntroduceDiracIndex::usage =
    "IntroduceDiracIndex[expr,k,list] locates all appearances of Dirac[___, \
\k,___] in expr and takes the specified momentum k out of the Dirac context \
to form the expression Dirac[___,Part[list,i],___] Scal[k,Part[list,i]]. This \
can be useful when one has to calculate Feynman integrals in a general way \
without being able to eliminate the Dirac structure by using DiracTrace. Note \
that the function only accepts indices that are declared as indices and that \
do not already appear in the given expression expr.";

 IntroduceScalIndex::usage =
    "IntroduceScalIndex[expr,k,list] locates all appearances of \
Scal[a_?MomentumQ, k_?MomentumQ] in expr and substitutes them with the \
expression Scal[a,Part[list,i]] Scal[k,Part[list,i]]. Note that squares of \
scalar products are not decomposed. Use IntroduceDoubleScalIndex for that \
matter. Scal[k,k] is not decomposed as well since one does not want that in \
\almost all the cases.";

 IntroduceDoubleScalIndex::usage=
    "IntroduceDoubleScalIndex[expr,k,list] locates one appearance of \
Scal[a_MomentumQ, k_?MomentumQ]^2 per summand in expr and substitutes them \
with the expression Scal[a,Part[list,1]] Scal[k,Part[list,1]] \
Scal[a,Part[list,2]] Scal[k,Part[list,2]]. The list must contain exactly two \
available indices (meaning indices that are declared as indices and that \
aren't already in expr.)";

 UsedIndices::usage =
    "UsedIndices[expr] gives a list of all indices (must be declared as \
indices with DeclareIndex) that appear in expr."
 
 
 DiracDecompose::usage = 
    "DiracDecompose[expr,variable] sorts expr according to all different \
Dirac[___]  in expr. It then saves the prefactors and the corresponding \
Dira[__] in variable[1,i] and variable[0,i] respectively (1<=i<=Number of \
different Dirac[___])."

 EpsilonSort::usage =
    "EpsilonSort[expr,reflist] sorts expressions in Epsilon[___] according
     to reflist.";

 EpsilonEpsilonContract::usage = 
    "EpsilonEpsilonContract[expr] contracts products of the form
     Epsilon[a,b,c,d] Epsilon[a,e,f,g]. Make sure that the same
     indices are the first in the List; use EpsilonSort otherwise.";

 OnShell::usage =
    "OnShell[{p1, m1}, {p2, m2}] returns a list with which the momenta in an
     expression can be set on shell: expr//.OnShell[{p1, m1}, {p2, m2}] does \
that.";

 SetOnShell::usage =
    "SetOnShell[{p1, m1}, {p2, m2}] sets Scal[p1,p1] = m1^2 etc.";

 ReplScal::usage = 
    "ReplScal[Scal[p1,p2], p1+p2+.. == ...] or similar generates a \
replacement
     list of the form {Scal[p1,p2] -> ....}. These rules are obtained by \
squaring 
     both sides of the equation. So 'p==p1+p2+q' will produce a different \
result
     than  '-q==p1+p2-p'; '-q-p1==p2-p' will return an empty list.";

 Sigma::usage =
    "Dirac[.. , Sigma[mu,nu], ...] can be used; Sigma[mu,nu] represents I/2 \
(Dirac[mu,nu] - Dirac[nu,mu]).
     Dirac[.. , Sigma[mu,nu], ...]/.SigmaRule inserts the explicit form."

  SigmaRule::usage =
    "Dirac[.. , Sigma[mu,nu], ...] can be used; Sigma[mu,nu] represents I/2 \
(Dirac[mu,nu] - Dirac[nu,mu]).
     Dirac[.. , Sigma[mu,nu], ...]/.SigmaRule inserts the explicit form."



 AddTraceRule::usage = "";

 ConditionalAddTraceRule::usage = "";

 ExtDiracTrace::usage = "";

 SortLst::usage = "";

 DiracPattern::usage = "";

 

(*********************************** Start the Package \
\*******************************************************)








Begin["`Private`"]






(*************************** Declare-Index, Mass, Momentum **** Queries   \
**********************)


(* Declaration of Momenta, masses and indices *)


ListOfIndices ={};  (* Keep an internal list of all the indices (what for?!) \
-> ContractAllIndices uses it *)

DeclareIndex[mu_] := (mu/: IndexQ[mu] = True ;
		      ListOfIndices= Union[ListOfIndices,{mu}] ;);
DeclareIndex[mu1_,mu2___] := ( DeclareIndex[mu1]; DeclareIndex[mu2]);
DeclareIndex[] :=ListOfIndices;

(* DeclareMass[m]: declares m as real as well *) 
DeclareMass[m_] := If[ m =!= 0 , m/: MassQ[m] = True; Conjugate[m] ^= m; ];
DeclareMass[m1_,m2___] := ( DeclareMass[m1] ; DeclareMass[m2]; );

DeclareMomentum[p_] := (p/: MomentumQ[p] = True; );
DeclareMomentum[p1_,p2___] :=( DeclareMomentum[p1] ; DeclareMomentum[p2]);



(* Queries for Momenta, ... *) 

MassQ[c_?NumberQ m_]:= MassQ[m]
MassQ[Plus[m1_,m2___]] := MassQ[m1] && MassQ[Plus[m2]];
NotMassQ[m_]:=!MassQ[m];


MomentumQ[c_?NumberQ p_]:= MomentumQ[p]
MomentumQ[Plus[p_,q___]] := MomentumQ[p] && MomentumQ[Plus[q]];
NotMomentumQ[p_]:=!MomentumQ[p];





       
DiracVectorQ[a_] := (IndexQ[a] || MomentumQ[a])
 (* L, R, Unit and Gamma5 are (treated as) Projectors here *)
    ProjectorQ[x_]:= Head[x]=== Projector;d





DeclarePolarizationVector[e_?AtomQ] := 
    (
     DeclareMomentum[e];
     PolarizationVectorQ[e] ^= True;
     Unprotect[MomentumQ];
     MomentumQ[Conjugate[e]] = True;
     Protect[MomentumQ];
     );


DeclarePolarizationVector[{e_,k_?MomentumQ}] := 
    (
     DeclareMomentum[e];
     PolarizationVectorQ[e] ^= True;
     Unprotect[MomentumQ];
     MomentumQ[Conjugate[e]] = True;
     Protect[MomentumQ];
     Scal[e,k] ^= 0;
     k /: Scal[Conjugate[e],k] = 0;
     );


DeclarePolarizationVector[e___,e1_] := 
    (
     DeclarePolarizationVector[e];
     DeclarePolarizationVector[e1];
     );





(******************* Useful stuff ********************)





OnShell[a_, b__] := Flatten[{OnShell[a], OnShell[b]}];
OnShell[{p_?MomentumQ, m_?MassQ}] = {Scal[p, p]  ->  m^2};
(* Allow m to be 0 [or any number]  *)
OnShell[{p_?MomentumQ, m_}]  = {Scal[p, p] -> m^2};



SetOnShell[e1_,e2___] := (SetOnShell[e1]; SetOnShell[e2];);
SetOnShell[{p_?MomentumQ,m_?MassQ}] := (p /: Scal[p,p] = m^2;);
(* Allow m to be 0 [or any number]  *)
SetOnShell[{p_?MomentumQ,m_}] := (p /: Scal[p,p] = m^2;);





ReplScal[Scal[p1_,p2_], leq_ == req_] := 
    Module[{eq, sol, sp},
	   eq = ( Scal[leq,leq]  ==  Scal[req,req] // DiracAlgebra);
	   eq = eq //. {Scal[p1,p2] -> sp};
	   sol = Solve[eq, sp][[1]];
	   sol = sol //. {sp -> Scal[p1,p2]};
	   Return[sol];
       ];




SigmaRule = {Dirac[a___,Sigma[mu_,nu_],b___] :> I/2 Dirac[a,mu,nu, b] -  I/2 \
Dirac[a,nu,mu, b]};




(******************************** One Rule for Dirac \
********************************)




Dirac[a___,0,b___] := 0;




(******************************* Scalar product \
**********************************)




SetAttributes[Scal,Orderless];
Format[Scal[p_,q_]]:=p\[CenterDot]q;
    Global`CenterDot = Scal;

    
    Scal[0,a_] := 0;
   
    
    ContractIndex[expr_,{mu_}]:= expr//.{Dirac[a___,mu,b___] Scal[c__,mu] :> \
Dirac[a,c,b],
                                         DiracTrace[Dirac[a___,mu,b___]] \
Scal[c__,mu] :> DiracTrace[Dirac[a,c,b]],
                                         Scal[a__,mu] Scal[c__,mu] :> \
Scal[a,c],
                                         Scal[a_?IndexQ,a_?IndexQ] :> d,
                                         Scal[mu,p_?DiracVectorQ]^2 -> \
Scal[p,p],
                                         (* conrtact indices in Dirac[] as \
well; this is a copy of a part of 'DiracAlgebra' *)
                                         Dirac[p1___,mu,mu,p2___] :> d \
Dirac[p1,p2],
                                         Dirac[p1___,mu,p2_,mu,p3___] :> \
(2-d) Dirac[p1,p2,p3],
                                         Dirac[p1___,mu,p2_,p3_,mu,p4___] :>
                                         4 Scal[p2,p3] Dirac[p1,p4]+(d-4) \
Dirac[p1,p2,p3,p4] /; 
                                           (IndexQ[p2] || MomentumQ[p2]) &&  \
(IndexQ[p3] || MomentumQ[p3]),
                                         (* If two indices are further apart: \
*)
                                         Dirac[p1___,mu,p2___,p3_,p4_,p5_,mu,\
p6___] :> 
                                         2 Dirac[p1,p5,p2,p3,p4,p6] - \
Dirac[p1,mu,p2,p3,p4,mu,p5,p6]
                                         };    
                                         (* this should not appear!! *)
                                         (* Power[Scal[a__,mu],n_?EvenQ] :> \
Power[Scal[a,a],n/2]  *)


    ContractIndex[expr_,{mu_,nu___}]:= \
ContractIndex[ContractIndex[expr,{mu}],{nu}];
    ContractIndex[expr_]:=expr;








(*******************************************
 Rules for Unit, L and R  
*********************************)





    L /: ProjectorQ[L] = True;
    R /: ProjectorQ[R] = True;
    Gamma5 /: ProjectorQ[Gamma5] = True;
    Unit /: ProjectorQ[Unit] = True;



    Unit /: Dirac[x___,Unit, y___] := Dirac[x,y];
    L/: Dirac[x___,L,L,y___]  := Dirac[x,L,y];
    L/: Dirac[x___,L ,R,y___] := 0;
    L/: Dirac[x___,a_?MomentumQ,L,y___] := Dirac[x,R,a,y] ;
    L/: Dirac[x___,a_?IndexQ,L,y___] := Dirac[x,R,a,y] ;
    L/: Dirac[x___,m_?MassQ,L,y___] := m Dirac[x,L,y] ;
    L/: Dirac[x___,Gamma5,L,y___] := - Dirac[x,L,y];
    L/: Dirac[x___,L,Gamma5,y___] := - Dirac[x,L,y];

    R/: Dirac[x___,R, R,y___] := Dirac[x,R,y];
    R/: Dirac[x___,R ,L,y___] := 0;
    R/: Dirac[x___,a_?MomentumQ,R,y___]:= Dirac[x,L ,a,y];
    R/: Dirac[x___,a_?IndexQ,R,y___]:= Dirac[x,L ,a,y];
    R/: Dirac[x___,m_?MassQ,R,y___]:= m Dirac[x,R,y];
    R/: Dirac[x___,Gamma5 , R,y___] := Dirac[x,R,y];
    R/: Dirac[x___,R , Gamma5,y___] := Dirac[x,R,y];

    Gamma5/: Dirac[Gamma5, Gamma5,a___] := Dirac[a];
    Gamma5/: Dirac[x___,a_?MomentumQ,Gamma5,y___] := - Dirac[x,Gamma5 ,a,y];
    Gamma5/: Dirac[x___,a_?IndexQ,Gamma5,y___]    := - Dirac[x,Gamma5 ,a,y]; 
    Gamma5/: Dirac[x___,a_?MassQ,Gamma5,y___]     :=   Dirac[x,Gamma5 ,a,y];

    (* And some rules for '-'  *)
    (* L/: *)    Dirac[a___,-L,b___] := - Dirac[a,L,b];
    (* R/: *)    Dirac[a___,-R,b___] := - Dirac[a,R,b];
    (* Gamma5/: *) Dirac[a___,-Gamma5,b___] := - Dirac[a,Gamma5,b];
   







    (********************** Linearity, Algebra and contractions \
********************************)



    (* Linearity for Dirac expressions *)

    DiracLinearity[expr_] := 
      Module[{res},
        (* Expand Everything within Dirac[] and Scal[] first *)
        res = expr //DiracScalExpand;
        res = res //.DiracLinearityList;
        Return[res]; 
      ];


    DiracLinearityList = 
      {  
        Dirac[x___,p_ + q_, y___] :> Dirac[x,p,y] + Dirac[x,q,y],
	Scal[x___,p_ + q_ , y___] :> Scal[x,p,y] + Scal[x,q,y],
        Dirac[a___, p_?MomentumQ  x_, b___]   :> x Dirac[a,p,b],
        Dirac[a___, p_?IndexQ  x_, b___]      :> x Dirac[a,p,b],
        Dirac[a___, p_?ProjectorQ  x_, b___]  :> x Dirac[a, p,b],
        (* Now the masses are not multiplied by Projectors or momenta *)
        Dirac[a___, m_?MassQ  x_, b___]       :> x Dirac[a,m,b],
	Dirac[a___, m_?MassQ , b___]          :> m Dirac[a,b],
        Scal[a_, p_?MomentumQ x_]             :> x Scal[a,p],
        Scal[a_, p_?IndexQ x_]                :> x Scal[a,p]
        (* Scal[a_, p_?MassQ x_]              :> p Scal[a,x] *)
      };






    (* Rules for Gamma expressions *)

    DiracAlgebraList = 
      {  
 
        Dirac[p1___,m_?MassQ,p2___]:> m Dirac[p1,p2],
        Dirac[p1___,Plus[m_?MassQ ,p2___] ,p3___] :> m Dirac[p1,p3] + \
Dirac[p1,Plus[p2],p3],
        Dirac[p1___,Times[e1___, m_?MassQ,e2___],p2___] :> m \
Dirac[p1,Times[e1,e2],p2] ,
   

        (* A few contractions for indices *)
        Dirac[p1___,mu_?IndexQ,mu_,p2___] :> d Dirac[p1,p2],
        Dirac[p1___,mu_?IndexQ,p2_,mu_,p3___] :> (2-d) Dirac[p1,p2,p3],
        Dirac[p1___,mu_?IndexQ,p2_,p3_,mu_,p4___] :>
          4 Scal[p2,p3] Dirac[p1,p4]+(d-4) Dirac[p1,p2,p3,p4] /; (IndexQ[p2] \
|| MomentumQ[p2]) &&  (IndexQ[p3] || MomentumQ[p3]),
        (* Dirac[p1___,mu_?IndexQ,p2_,p3_,p4_,mu_,
          p5___] :> -2 Dirac[p1,p4,p3,p2,p5]+(4-d) Dirac[p1,p2,p3,p4,p5], *)
        (* If two indices are further apart: *)
        Dirac[p1___,mu_?IndexQ,p2___,p3_,p4_,p5_,mu_,p6___] :> 
          2 Dirac[p1,p5,p2,p3,p4,p6] - Dirac[p1,mu,p2,p3,p4,mu,p5,p6],

        (* The same contractions for momenta *)
        Dirac[p1___,p_?MomentumQ,p_,p2___] :> Scal[p,p] Dirac[p1,p2] ,
        Dirac[p1___,p_?MomentumQ,p2_,p_,p3___] :> 2 Scal[p,p2] Dirac[p1,p,p3] \
- Scal[p,p] Dirac[p1,p2,p3] /; (IndexQ[p2] || MomentumQ[p2]),
        (* If two momenta are further apart: *)
        Dirac[p1___,p_?MomentumQ,p2___,p3_,p_,p4___] :> 2 Scal[p,p3] \
Dirac[p1,p,p2,p4] - Dirac[p1,p,p2,p,p3,p4]/; (IndexQ[p3] || MomentumQ[p3])

      };





    (* DiracAlgebra itself *)

    DiracAlgebra[expr_] := 
      Module[{res},
        (* Expand Everything within Dirac[] and Scal[] first *)
        res = expr //DiracScalExpand;
        res = res // DiracLinearity;
        res = res //.DiracAlgebraList;
        Return[res]; 
      ];


   



    (* Contractions for all indices  *)

    ContractAllIndices[expr_] :=  Module[{tmp=expr},
      tmp = tmp//DiracAlgebra//Expand;
      For[i=0,i<=Length[ListOfIndices],i++,
          tmp=ContractIndex[tmp,ListOfIndices]];
          Return[tmp]];

ListOfIndices2={lor1,lor2,lor3,lor4,lor5,lor6,lor7,lor8};


 ContractAllIndices2[expr_] :=  Module[{tmp=expr},
      tmp = tmp//DiracAlgebra//Expand;
      For[i=0,i<=Length[ListOfIndices2],i++,
          tmp=ContractIndex[tmp,ListOfIndices]];
          Return[tmp]];





    (************************************** Cosmetic functions \
**************************************)




    DiracSort[expr_,RefList_]:=Module[{res,RL},
		RL = Join[{R},{L},{Gamma5},RefList];
		(* These are on the left anyway *)
		res= expr;
		res = res//.{Dirac[a___,b_?DiracVectorQ,c_?DiracVectorQ,dd___] :> - \
Dirac[a,c,b,dd] + 2 Scal[b,c] Dirac[a,dd] /; 
                      (Position[RL,b][[1,1]] > Position[RL,c][[1,1]]  )};
		Return[res]
               ];
               (* What happens for Dirac[___,Conjugate[eps1],___,eps1,___]? \
Check note for EpsilonSort *)


  
    (* A routine to Collect and simplify Dirac[] expressions *)
    DiracCollect[expr_] := Collect[expr,Dirac[___]];
    DiracCollect[expr_,func_] := Collect[expr,Dirac[___],func];



    (* one alias *)
    DiracFactor[expr_] := DiracCollect[expr,Factor];



    (* Expand everything inside Dirac and Scal *)
    DiracScalExpand[expr_] := expr//.{Dirac[a___]:> Expand/@Dirac[a], \
Scal[a___]:> Expand/@Scal[a] };






    (************************************** Dirac equation and Reference List \
*******************************)


    (* Rules to use the Dirac equation *)
    (* Instead of the distinction U and V could simply use m for particles \
and -m for antiparticles  *)


    (*
       UBarLeft::usage =
         "UBarLeft[expr,p,m] turns an expression of the form Dirac[p,...]
          into m Dirac[...] as if the Dirac eq had been used.";

       VBarLeft::usage =
         "VBarLeft[expr,p,m] turns an expression of the form Dirac[p,...]
          into -m Dirac[...] as if the Dirac eq had been used.";

       URight::usage =
         "URight[expr,p,m] turns an expression of the form Dirac[...,p]
          into m Dirac[...] as if the Dirac eq had been used.";

       VRight::usage =
         "VRight[expr,p,m] turns an expression of the form Dirac[...,p]
          into -m Dirac[...] as if the Dirac eq had been used.";
     *)




    URight[expr_,p_,m_] := expr//.{Dirac[x___,p]-> m Dirac[x]};
    VRight[expr_,p_,m_] := expr//.{Dirac[x___,p]-> - m Dirac[x]};
 
    UBarLeft[expr_,p_,m_] := expr//.{Dirac[p,x___]-> m Dirac[x],
				Dirac[L,p,x___]->    m Dirac[R,x],
				Dirac[R,p,x___]->    m Dirac[L,x],
				Dirac[Gamma5,p,x___]-> - m Dirac[Gamma5,x]};
    VBarLeft[expr_,p_,m_] := expr//.{Dirac[p,x___]-> - m Dirac[x],
				Dirac[L,p,x___]->    - m Dirac[R,x],
				Dirac[R,p,x___]->    - m Dirac[L,x],
				Dirac[Gamma5,p,x___]->   m Dirac[Gamma5,x]};
  



    (* Produce a Reference list to order expressions in Dirac[];
       Give leftmost and rightmost Momentum or Index;
       UseDiracEquation needs this function.                    *)
    MakeReflist[expr_,{left_,right_}] :=
      Module[{mom,in,lst},
	  lst=Cases[expr,Dirac[___],{0,Infinity}];
	  lst=Flatten[lst/.Dirac->List];
	  lst = Union[lst];
	  (* Get rid of R and L and whatever *)
	  mom=Select[lst,MomentumQ ];
	  in = Select[lst,IndexQ ];
	  lst = Union[in,mom];
	  (* Drop the elements you need *)
	  lst = DeleteCases[lst,left];
	  lst = DeleteCases[lst,right];
	  (* And add them in the right place *)
	  lst = Append[lst,right];
	  lst = Prepend[lst,left];
	  Return[lst];
           ];

    (* Produce a Reference list to order expressions in Dirac[];
       Give leftmost  Momentum or Index;
       UseDiracEquation needs this function.                    *)
    MakeLeftReflist[expr_,{left_}] :=
      Module[{mom,in,lst},
	  lst=Cases[expr,Dirac[___],{0,Infinity}];
	  lst=Flatten[lst/.Dirac->List];
	  lst = Union[lst];
	  (* Get rid of R and L and whatever *)
	  mom=Select[lst,MomentumQ ];
	  in = Select[lst,IndexQ ];
	  lst = Union[in,mom];
	  (* Drop the elements you need *)
	  lst = DeleteCases[lst,left];
	  (* And add them in the right place *)
	  lst = Prepend[lst,left];
	  Return[lst];
           ];


    (* Produce a Reference list to order expressions in Dirac[];
       Give rightmost and rightmost Momentum or Index;
       UseDiracEquation needs this function.                    *)
    MakeRightReflist[expr_,{right_}] :=
      Module[{mom,in,lst},
	  lst=Cases[expr,Dirac[___],{0,Infinity}];
	  lst=Flatten[lst/.Dirac->List];
	  lst = Union[lst];
	  (* Get rid of R and L and whatever *)
	  mom=Select[lst,MomentumQ ];
	  in = Select[lst,IndexQ ];
	  lst = Union[in,mom];
	  (* Drop the elements you need *)
	  lst = DeleteCases[lst,right];
	  (* And add them in the right place *)
	  lst = Append[lst,right];
	  Return[lst];
           ];



    (* Sort and use Dirac equation *)
    UseDiracEquation[{LeftMom_,LeftMass_},expr_,{RightMom_,RightMass_}] := 
      Module[{res,RefList},
	  (* Need to build a reflist first *)
	  RefList=MakeReflist[expr,{LeftMom,RightMom}];
	  res = DiracSort[expr,RefList];
	  res = DiracCollect[res];
	  res = UBarLeft[res,LeftMom,LeftMass];
	  res = URight[res,RightMom,RightMass];
	  res = DiracCollect[res];
	  Return[res];
           ];


    UseDiracEquation[{},expr_,{RightMom_,RightMass_}] := 
      Module[{res,RefList},
	  (* Need to build a reflist first *)
	  RefList=MakeRightReflist[expr,{RightMom}];
	  res = DiracSort[expr,RefList];
	  res = DiracCollect[res];
	  res = URight[res,RightMom,RightMass];
	  res = DiracCollect[res];
	  Return[res];
           ];


    UseDiracEquation[{LeftMom_,LeftMass_},expr_,{}] := 
      Module[{res,RefList},
	  (* Need to build a reflist first *)
	  RefList=MakeLeftReflist[expr,{LeftMom}];
	  res = DiracSort[expr,RefList];
	  res = DiracCollect[res];
	  res = URight[res,LeftMom,LeftMass];
	  res = DiracCollect[res];
	  Return[res];
           ];



    (***************************************** Dirac Adjunction \
*****************************************)




    DiracAdjunction[L] := R;
    DiracAdjunction[R] := L;
    DiracAdjunction[Gamma5] := -Gamma5;
    DiracAdjunction[Unit] := Unit;


    DiracAdjunction[a_ b_] := 
      DiracAdjunction[a] DiracAdjunction[b];

    DiracAdjunction[a_ + b_] := DiracAdjunction[a] + DiracAdjunction[b];



    (* Indices, momenta and masses remain the same (p represents p-slash; mu \
represents gamma_mu, so ... )
       but not polarizationvectors! *)
    DiracAdjunction[a_] := Conjugate[a] /; PolarizationVectorQ[a]; 
    DiracAdjunction[a_] := a /; (DiracVectorQ[a] || MassQ[a]) && \
PolarizationVectorQ[a] =!= True; 
 
    DiracAdjunction[a_] := Conjugate[a] /; FreeQ[a,Dirac] && FreeQ[a,Scal];

    DiracAdjunction[Dirac[a___]] := 
      Dirac[(DiracAdjunction/@Reverse[{a}])/.{List -> Sequence}];

    DiracAdjunction[Scal[a_,b_]] := Scal[DiracAdjunction[a], \
DiracAdjunction[b]];

    (* Can also have Scal[a,b]^n *)

    DiracAdjunction[Scal[a_,b_]^n_Integer] := Scal[DiracAdjunction[a], \
DiracAdjunction[b]]^n;

    (* And 1/(...) as well ... *)

    DiracAdjunction[A_^n_Integer] := DiracAdjunction[A]^n;


    (********************************************  Traces \
***********************************************)




    (* A Routine for Gamma Traces *)


    (* Need Epsilon (for Gamma5)  *)

    (* Convention \epsilon^{0123} = -1 *)


    Epsilon[a___, 0, b___] := 0;
    Epsilon[a___,b_,c___,b_,e___] := 0;


    (* Contract and Sort Epsilons *)

    ContractScalEps[expr_] := 
      expr//.{Scal[a_?IndexQ,b_] Epsilon[c___,a_,dd___] :> Epsilon[c,b,dd]};

    EpsilonSort[expr_,reflist_] := 
        Module[{res},
	       res=expr//.{Epsilon[a___,b_,c_,dd___] :> - Epsilon[a,c,b,dd] /; 
			   ( Position[reflist,b,{1}][[1,1]] > Position[reflist,c,{1}][[1,1]])};
	       Return[res]];
    (* Need the level-specification in Position, otherwise e.g. reflist = \
{Conjugate[e1],e1} goes wrong! *)


    EpsilonContract1 = { Epsilon[mu_, nu_, rh_, si_] Epsilon[mu_, be_, ga_, \
de_]  :> 
       - Scal[be,nu] Scal[ga,rh] Scal[de,si] - Scal[ga,nu] Scal[de,rh] \
Scal[be,si] 
       - Scal[de,nu] Scal[be,rh] Scal[ga,si]
       + Scal[be,nu] Scal[de,rh] Scal[ga,si] + Scal[ga,nu] Scal[be,rh] \
Scal[de,si]
       + Scal[de,nu] Scal[ga,rh] Scal[be,si] };


    EpsilonContract2 = { Epsilon[mu_, nu_, rh_, si_] Epsilon[mu_, nu_, ga_, \
de_]  :> 
       - 2 ( Scal[ga,rh] Scal[de,si] - Scal[de,rh] Scal[ga,si] ) };

    EpsilonContract3 = { Epsilon[mu_, nu_, rh_, si_] Epsilon[mu_, nu_, rh_, \
de_]  :> 
       - 6 Scal[de,si] };

    EpsilonContract4 = { Epsilon[mu_, nu_, rh_, si_] Epsilon[mu_, nu_, rh_, \
si_]  :>
	- 24 };

    EpsilonEpsilonContract[expr_] := 
	Module[{res},
		res = expr//.EpsilonContract4;
		res = res//.EpsilonContract3;
		res = res//.EpsilonContract2;
		res = res//.EpsilonContract1;
		Return[res];
	  ];



    (* DiracTraceList can handle expressions up to four Gammas only *)

    TraceMoreThanFour := 
      { (* Make sure there's no projector in front of the expression! *)
        DiracTrace[Dirac[p_?DiracVectorQ,ex__]] :> 
          Module[{res, expr, i, k, DiracTmp}, expr = DiracTmp[p,ex];
          res = Sum[(-1)^k Scal[expr[[1]], expr[[k]]] *
                DiracTrace[ Delete[expr, {{1}, {k}}] ], {k, 2, \
Length[expr]}];
          res = res //. {DiracTmp -> Dirac};
          res]}


    TraceMoreThanFourGamma5 := 
      {
        DiracTrace[Dirac[Gamma5, mu_, a__]] :> 
          Module[{res},
                 DeclareIndex[kk1, kk2, kk3];
                 res = -I/6  Epsilon[mu, kk1, kk2, kk3]  DiracTrace[ \
Dirac[kk1, kk2, kk3, a]]; 
                 res]
      };

    TraceGammaRL := 
      {
        DiracTrace[Dirac[R,a___]] :> 1/2 DiracTrace[Dirac[a]] + 1/2 \
DiracTrace[Dirac[Gamma5,a]],
        DiracTrace[Dirac[L,a___]] :> 1/2 DiracTrace[Dirac[a]] - 1/2 \
DiracTrace[Dirac[Gamma5,a]]
      };


    (* Rules for Traces up to four Gammas *)

    DiracTraceList = 
        { 
          DiracTrace[Dirac[b_,a___]] :> 0 /; EvenQ[ Length[{a}]] && \
Not[ProjectorQ[b]],
	  DiracTrace[Dirac[b_?ProjectorQ,a___]] :> 0 /; OddQ[ Length[{a}]],
          DiracTrace[Dirac[]] -> 4,
	  DiracTrace[Dirac[mu_?DiracVectorQ, nu_?DiracVectorQ]] :> 4 Scal[mu, nu], 
	  DiracTrace[Dirac[mu_?DiracVectorQ, nu_?DiracVectorQ, rh_?DiracVectorQ, 
	      si_?DiracVectorQ]] :> 
		  4 (Scal[mu, nu]Scal[rh, si] - 
		     Scal[mu, rh]Scal[nu, si] + Scal[mu, si]Scal[nu, rh]),
          DiracTrace[Dirac[Gamma5]] :> 0,
          DiracTrace[Dirac[Gamma5,mu_?DiracVectorQ, nu_?DiracVectorQ]] :> 0,
          DiracTrace[Dirac[Gamma5,mu_?DiracVectorQ, nu_?DiracVectorQ, \
rh_?DiracVectorQ, 
	      si_?DiracVectorQ]] :> 4 I Epsilon[mu, nu, rh, si],
          DiracTrace[Dirac[R]] :> 2,
          DiracTrace[Dirac[R,mu_?DiracVectorQ, nu_?DiracVectorQ]] :> 2 \
Scal[mu,nu],
          DiracTrace[Dirac[L]] :> 2,
          DiracTrace[Dirac[L,mu_?DiracVectorQ, nu_?DiracVectorQ]] :> 2 \
Scal[mu,nu]
		  };

    (* This can be filled at runtime with more rules to make it faster *)
    (* Make sure that DiracTraceLinearity rules are already applied to this \
list!! *)

    ExtendedDiracTraceList = {};



    DiracTraceLinearity[expr_] := Module[{res}, 
      res = expr//.{ DiracTrace[a___]:> Expand/@DiracTrace[a]};
      res = res//.{DiracTrace[Plus[a_,b__]] :> DiracTrace[a] + \
DiracTrace[Plus[b]],
      DiracTrace[a_ Dirac[b___]] :> a DiracTrace[Dirac[b]]}; 
      Return[res];];
    



    (* Here things like DiracTrace[Dirac[...]] enter *)
    DiracTraceAlgebra[expr_] :=
        Module[{res},
           res = expr;
           res = res//.SigmaRule;
           res = res//DiracAlgebra;
	   res = res//DiracTraceLinearity;
	   res = res//.ExtendedDiracTraceList;
           (* res = res//DiracTraceLinearity; *)
	   res = res//.DiracTraceList;
           (* res = res//DiracTraceLinearity; *)
           (* For everything that's not traced yet; missing: L,R, Gamma5 *)
           While[ Not[FreeQ[res,DiracTrace,{0,Infinity}]],
             res = res//.TraceGammaRL;
             res = res//.TraceMoreThanFourGamma5;
	     res = res//.TraceMoreThanFour;
	     res = res//.ExtendedDiracTraceList;
	     res = res//.DiracTraceList;
	     res = Expand[res];
	     res = res//ContractScalEps;
                ];
	       Return[res];
	   ];

   (* How to add rules into ExtendedDiracTraceList *)



    DiracPattern[expr_] :=Module[{res,P,PT},
       	res= expr/.{ a_?DiracVectorQ :> PT[P[a,Blank[]], DiracVectorQ]};
       	res = res //.{P -> Pattern, PT -> PatternTest };
       	Return[res];];


    SortLst[expr_]:= Module[{res},
	res = expr/.{Dirac->List,DiracTrace -> List} // Flatten;
	res =Complement[res,{R,L,Gamma5}];
	Return[res];
		];


    ExtDiracTrace[DiracTrace[Dirac[expr___]]] := Module[{res},
	res =  DiracTraceAlgebra[DiracTrace[Dirac[expr]]];
	res = EpsilonSort[res,SortLst[{expr}]];
	res  = Factor[res];
	Return[res];];


    AddTraceRule[expr_] :=  Module[{tmp,pat},
	tmp = ExtDiracTrace[expr];
        pat = DiracPattern[expr];
        ExtendedDiracTraceList = Union[ExtendedDiracTraceList , {pat -> \
tmp}];
		];


    (* Adds a rule only, if there is no rule yet *)
    ConditionalAddTraceRule[expr_] := Module[{},
	If[ DiracTraceRuleExists[expr] , Goto[exit] ];
	AddTraceRule[expr];
	Label[exit];
	];



    DiracTraceRuleExists[DiracTrace[Dirac[a___]]] := 
      Module[{},
             (* If Dirac[] is entered *)
             If[a == Sequence[], Return[True]];
             (* These are zero anyway *)
             If[OddQ[Length[Dirac[a]]] && Not[ProjectorQ[Dirac[a][[1]] ] ], \
Return[True] ];
             If[EvenQ[Length[Dirac[a]]] && ProjectorQ[Dirac[a][[1]]], \
Return[True] ];
             (* Check if there's a rule *)
             If[(DiracTrace[Dirac[a]]//.DiracTraceList) =!= \
DiracTrace[Dirac[a]] , Return[True] ];
	     If[(DiracTrace[Dirac[a]]//.ExtendedDiracTraceList) =!= \
DiracTrace[Dirac[a]] , Return[True] ];
             Return[False];
            ];
        


    DiracExprMinusTwo[DiracTrace[Dirac[a___]]] := 
        Module[{},
                If[a == Sequence[], Return[DiracTrace[Dirac[]]]; ];
                If[Length[{a}] <= 2 , Return[DiracTrace[Dirac[a]]]; ];
                Return[DiracTrace[Dirac[Drop[{a}, -2]] /. {List -> \
Sequence}]];
              ];


    (* Recursively learn new rules; What happens for Dirac[___, \
Conjugate[e1], b___] ?! *)
    LearnDiracTraceRule[DiracTrace[Dirac[a___]]] := 
        Module[{GEMT},
                  If[DiracTraceRuleExists[DiracTrace[Dirac[a]]], Return[]; ];
                  GEMT = DiracExprMinusTwo[DiracTrace[Dirac[a]]];
                  If[DiracTraceRuleExists[GEMT], \
ConditionalAddTraceRule[DiracTrace[Dirac[a]]]; 
                         Return[]; ];
                  LearnDiracTraceRule[GEMT];
                  ConditionalAddTraceRule[DiracTrace[Dirac[a]]];
              ];











    (******************************* Squaring Dirac-Expressions \
*************************************)






    (* Helper function: DiracProduct[a Dirac[b], ___]  Joins all the Dirac[] \
expressions *)

    (* 'Linearity' *)
    DiracProduct[a___, b_ + c_ , dd___] := DiracProduct[a,b,dd] + \
DiracProduct[a,c,dd];
    DiracProduct[c___, a_ Dirac[b___],dd___] := a \
DiracProduct[c,Dirac[b],dd];


    (* Expand argument, if there is a Dirac in there; that is if there is \
more than a Dirac in there... *)
    DiracProduct[c___,a_,dd___]:= 
      DiracProduct[c,Expand[a],dd] /; (Not[FreeQ[a,Dirac]] && Not[Head[a] === \
Dirac]);

    (* Join *)
    DiracProduct[Dirac[a___],Dirac[b___],c___] := DiracProduct[Dirac[a,b],c];

    (* Get rid of the 'DiracProduct' if the work is done *)
    DiracProduct[Dirac[a___]] := Dirac[a];

    (* And: *)
    DiracProduct[a___,0,b___] := 0;





    (* DiracSquare itself: *)
    (* DiracSquare[expr,one,two] returns the trace of expr * one * \
DiracAdjunction[expr] * two *)

    DiracSquare[expr_,one_,two_] := Module[{res},
					 res= DiracProduct[expr,one,DiracAdjunction[expr],two];
					 res = res//DiracAlgebra;
                                         res = res//DiracTrace;
					 res = DiracTraceAlgebra[res];
					 Return[res];
				     ];


    (******************************* Added functions by Kay Bieri \
**********************************)

UsedIndices[expr_]:=
Module[{usedlist},
      usedlist={};
      For[i=1,i<Length[ListOfIndices]+1,i++,
        usedlist=
          Append[usedlist,
            Cases[expr,Part[ListOfIndices,i],{0,\[Infinity]}]]];
      usedlist=Intersection[Evaluate[Flatten[usedlist]],ListOfIndices];
      Return[usedlist];
      ];

IntroduceDiracIndex[abra_,momentum_?MomentumQ,mylist_?ListQ]:=
    Module[{emacs,newlist={},removelist={},res=abra,a,b,c},
      If[Length[ListOfIndices]==0,
          Return[res],
          removelist=Complement[ListOfIndices,UsedIndices[res]];
          newlist=Intersection[Evaluate[Flatten[mylist],removelist]];
          If[Length[newlist]==0,
            Return[res],
            res=res*emacs//Expand;
            For[i=1,i<Length[newlist]+1,i++,
              res=res/.{Dirac[a___,momentum,b___]*c_->
                      c*Scal[Part[newlist,i],momentum] Dirac[a,
                          Part[newlist,i],b]}];
            emacs=1 ;
            Return[res]]];
      ];
IntroduceDiracIndex[abra_,momentum_,list_]:=abra;


IntroduceScalIndex[abra_,momentum_?MomentumQ,mylist_?ListQ]:=
    Module[{saved,emacs,newlist={},removelist={},res=abra,a,b,c},
      If[Length[ListOfIndices]==0,
          Return[res],
          removelist=Complement[ListOfIndices,UsedIndices[res]];
          newlist=Intersection[Evaluate[Flatten[mylist],removelist]];
          If[Length[newlist]==0,
            Return[res],
            res=Evaluate[res //. Scal[momentum,momentum] -> saved];
            res=res*emacs//Expand;
            For[i=1,i<Length[newlist]+1,i++,
              res=res/.{Scal[a_?MomentumQ,momentum]*c_->
                      c*Scal[Part[newlist,i],momentum] Scal[a,
                          Part[newlist,i]]}];
            emacs=1;
            res=res/.saved -> Scal[momentum,momentum];
            Return[res]]];];
IntroduceScalIndex[abra_,momentum_,list_]:=abra;



IntroduceDoubleScalIndex[abra_,momentum_?MomentumQ,mylist_?ListQ]:=
    Module[{emacs,squ,newlist={},removelist={},res=abra,a,b,c},
    Clear[emacs,squ,a,b];
    If[Length[ListOfIndices]==0,Print["List of indices is empty"];
          Return[res],
          removelist=Complement[ListOfIndices,UsedIndices[res]];
          newlist=Intersection[Evaluate[Flatten[mylist],removelist]];
          If[Length[newlist]==2,
            res=res//.{Scal[momentum,momentum]^2->squ};
            res=res*emacs//Expand;
            
            res=res/.{Scal[a_?MomentumQ,momentum]^2*c_->
                    Scal[Part[newlist,1],momentum] Scal[a,Part[newlist,1]]*
                      Scal[Part[newlist,2],momentum]*Scal[Part[newlist,2],a]*
                      c};
            res=res/.squ->Scal[momentum,momentum]^2;
            emacs=1,
            Print["There are not exactly two allowed indices specified!"]];
          Return[res]];];

IntroduceDoubleScalIndex[abra_,momentum_,list_]:=abra;

DiracDecompose[expr_,variable_]:=
    Module[{i=0,len=0,zres=expr,mist=0},
    Clear[variable];
      mist=Cases[expr,Dirac[___],{0,\[Infinity]}]//Union;
      len=Length[mist];
      If[len==0,Return[zres],zwischen=Collect[zres,mist];
        While[i<len,i++;
          variable[0,i]=Part[mist,i];
          variable[1,i]=Coefficient[zwischen,Part[mist,i]]];
        Print["Number of different Dirac[___] - terms: "];
        Return[len];
        ];
      ];


    (*************************************************************************\
**********************)






    Protect["Fermions`*"];


    (* Things that have to be writable at runtime *)

    Unprotect[ListOfIndices,ExtendedDiracTraceList];



  End[];    (* "`Private`" *) 

EndPackage[];



