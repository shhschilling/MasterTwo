(* ::Package:: *)

BeginPackage["Integrals`",{"Fermions`"}]
Unprotect["Integrals`*"];


IntegralsInfo[] := (
		    Print[" \n Integrals is a package that can perform two dimensional Tensor reduction for up to 4 free lorentzindices, it can perform onedimensional Tensor reduction for up to 8 free indicices. It can integrate integrals with integrands of the form 1/(q1^2-m1^2)^a(q2^2-m2^2)^b(q1+q2)^c."];
		    Print["\n Available commands:"];
		    Print[" (For any of these commands '?command' should help.)\n"];
		    Print[Names["Integrals`*"]]);


     
     (* Note also that we use the convention 'd = 4 -2eps', where d is the space-time dimension. *)
    
     

(************************    USAGE                                      ***************************)
(******************Headers of MasterTwo*********************************)
        AD::usage="AD is the header of a list of  propagators. For two-loop integrals it is the header of factorising one-loop integrals after the usage of SimplifyPropagator."
       den::usage="is the header of a propagator. den[q1,m1] corresponds ...
           to  1/(Scal[q1,q1]-m1^2)."
       
   G::usage="G is the header of a list of  propagators for two-loop ...
            integrals with one vanishing mass. G[i[M1,q1],i[M2,a2],i[0,a3]]  corresponds to  
     1/((q1^-M1^2)^a1(q2^-M2^2)^a2(q1+q2)^a3) with at least one mass in the propagators 0 and q1 and q2 being loop momenta being declared with DeclareLooopMomentum."
       
       
       
        
       i::usage="i is the header of propagators in a list of propagators of the form  G[i[M1,q1],i[M2,a2],i[0,a3]] corresponding to  1/(q1^-M1^2)^a1(q2^-M2^2)^a2((q1+q2))^a3] with at least one mass in the propagators 0 and q1 and q2 being loop momenta being declared with DeclareLooopMomentum."
 Ne::usage="Header of prefactors of the form Ne[m]=(mu^2 /m ^2)^eps 2 ^{2 eps}\pi^ eps \Gamma(1+eps) as defined in  hep-ph/9910220."
       N2::usage="Header of prefactors of the form N2[m]=Ne[m]^2=(mu^{2} /m ^{2})^{2eps} 2 ^{4 eps}\pi^{2 eps} \ Gamma(1+eps)^2 as defined in  hep-ph/9910220."

       log::usage="Abbreviation for natural Logarithms Log to prevent Magthematica to simplify expressions containg log"       
       PoLi2::usaage="PoLi2[x] is the header of the polilogarithm ...
           Li_2(x)=\integral_0^1 \f {\ln (t)} {t-1/x}dt. Its definition corresponds to the Mathematica-Function PolyLog[2,x] of the function x."
       
      
       
       ClausenCl2::usage =
      "ClausenCl2[x] corresponds to the Mathematica-function
          ClausenCl[2,x]. Introduced to prevent Mathematica from ...
          further simplications. "
       
       ClausenCl::usage =
      "ClausenCl[n,x] gives the Clausen function of order n."

       
       SUNT::usage="SUNT is the header of the colour generator of SU(3)_c,SUNT[a,i,j] stands for T^a_{ij}SUNT[a,b] stands
           for the product of the colour generators T^a T^b. "
       SUNF::usage="SUNF is the header of the structure. SUNF[a,b,c] stands for f_{a b c}." 

       
       
       
       
       scalstructure::usage=""
       
       scalerule::usage=""
       fc::usage=""
       denrule::usage=""
       denruleback::usage=""
       xrule::usage=""
       colorrules::usage=""
       deletprop::usage=""
       Deletprop::usage=""
       delet::usage=""
       forwardpower::usage=""
       middlepower::usage=""
       backpower::usage=""
       power::usage=""
       subrule::usage=""
       FinalPart::usage=""
       Denrule::usage=""
       DenruleOneBack::usage=""
       DenruleOne::usage=""
       SimplifyPropagator::usage="Brings propagators of integrands to the form needed for loop integration. Example:
AD[den[q1,m1], den [q1,m1], den[q1,m1], den[q2,m2], den[q2,m2], den[q3,m3] is transformed to G[m1,3], den [m2,2], [m3,1]."
       DenruleBack::usage=""
       
	 
           
       Color::usage="replaces color structures depending on generators and structure constants of SU(3)_c on expressions only depending on generators or scalars."
     

       Translation::usage=""
FeynArtsToMasterTwo::usage="Translates the output of FeynArts into the notation needed for MasterTwo."





       TwoZero::usage=""

       TwoEqual::usage=""   

       

      

       ScalIntOne::usage="ScalIntOne[AD[i[m,n]]]allows the calculation
           of scalar one loop integrals by replacing the propagator
           structure AD[i[m,n]]  with (pi^2/((m^2)^ {n-2})) Ne(m)
           C^{(1)}_n. Not that in order to obtain the standard 
           Integral (\mu ^{2eps} /({2 Pi}^{2D})) {q^2-m^2}^{-n}  
           the result hast to be multiplied by 1/((2 Pi)^4)."

ScalIntTwo::usage= "ScalIntTwo[G[i[M1, a1], i[M2, a2], i[0, a3]
calculatesscalar twoloop integrals of the form  nu^(4 eps)/(2 Pi^(-4 eps)) Integrate[( d^q1
 d^q_2)/((2 Pi)^(2D)) *1/(q1^-M1^2)^a1(q2^-M2^2)^a2((q1+q2))^a3]
 with at least one mass in the propagators 0 by replacing the
propagator structure G[i[M1, a1], i[M2, a2], i[0, a3]] with the analytical result of the  
corresponding integral.A prefactor N2[M1] as defined in  
hep-ph/9910220, equation (61) is introduced. Results are Taylor  
expanded in epsilon until 0. order. Loop momenta  
have to be declared with DeclareLoopMomentum. Before usage the  
integrands have to be transformed via substitutions in such a  
way that the same masses appear always with the same loop momenta  
with the help of the function SimplifyPropagator in all arising  
diagramms.  A factor x1=M2^2/M1^2 is introduced to shorten the result. Prefactors  
are expressed proportional to N2[M1] defined in eq. (65) of the  
manual. The simplifications of the prefactors is performed.
ScalIntTwo can also handle factorizing integrals of the form AD[i[M1, a1],i[M2,a2]]."


ScalIntTwoFull::usage= "ScalIntTwoFull[G[i[M1, a1], i[M2, a2], i[0, a3]]]: Functions as ScalIntTwo, but sorts the denominators before integration.A flag \.b4model\.b4 for the chosen model (SM or 2hdm) can be set. If model=SM,  MW is replaced by M1, if model=2hdm, MH is replaced by M1."


ScalIntTwoThreeMasses::usage = "ScalIntTwoThreeMasses[G[i[M1, a1], i[M2, a2], i[M1, a3]]] allows to
calculate scalar loop integrals independend of external momenta
and with up two different masses by  replacing propagator structures
G[i[m1, n1], i[m2, n2], i[m1, n3]] with the analytical result for 
scalar twoloop integrals as given in eq. (63) of the manual. 
Note that you have to give an explicit value x>0 for the mass
relation betwen x=M2^2/M1^2. The final 
results is Taylor expanded in eps up to second order."





facruleone::usage = "calculates standard one-loop integrals
AD[den[q,m],...,den[q,m]] (a propagator terms) by making the replacement  I*(((-1)^a*Pochhammer[1 + eps, a - 3]*Ne[m1]*
       (m1^2)^(2 - a))/(a - 1)!).. Note that in order to obtain standard one loop integrals
\int d^dp \nu^{2 eps } /(2 Pi)^D 1/(q^2-m^2)^a the final result has to be multiplied by 1/(16 Pi^4)."


subden::usage = "is a replacement rule which reorders propagators of the form
AD[l___,den[a_,m_],r__] in the form needed for the loop integration."



GammaExpand::usage= "GammaExpand[a] return Gamma[Expand[a]]."

GammaSimplify::usage = "GammaSimplify[expr] checks the list GammaSimplifyList and returns expr //. GammaSimplifyList. Rules to this list can be added in the package."

BetaToGamma::usage = "BetaToGamma searches an expression for Beta-functions and substitutes them with Gamma-functions."

ReplOneMinus::usage = " Expr //. ReplOneMinus[u] replaces 1-u by OneMinus[u] and -1+u by -OneMinus[u]. Afterwards it checks whether there are expressions like (OneMinus[u] u)^q. These expressions are then broken apart to OneMinus[u]^q u^q. Note that instead of u you can have a list of parameters like ReplOneMinus[{u,v,...}] or just ReplOneMinus[u,v,...]. Then the replacing rules are applied to all given variables."



BetaIntegrals::usage = " Expr //. BetaIntegrals[u] replaces  'u^p_.  OneMinus[u]^q_.' by 'Beta[p + 1, q + 1]'. You can do this simultaneously for several variables by using BetaIntegrals[{u,v,...}] or just BetaIntegrals[u,v,...]. Make sure to always use '//.' (ReplaceRepeated) and not '/.' (ReplaceOnce)!"

BetaRemains::usage = "BetaRemains[u] works just like BetaIntegrals and substitutes 'u^p_.  OneMinus[u]^q_.' with 'Beta[p + 1, q + 1]'. The difference is that it does so only for either p or q equal to zero (therefore the name '-Remains'.). Be sure to always use 'BetaIntegrals' first. You might get false expressions otherwise."

RemoveOddPowersOfIntVar::usage = "RemoveOddPowersOfIntVar[expr,var] removes all odd powers of var in expr (assuming that expr is the numerator of the Feynman integral about to be calculated). Make sure that the denominator of the Feynman ntegral is of the form (var^2-C)^a."

TenScal::usage = "";

EightScal::usage = "";

SixScal::usage = "";

FourScal::usage = "";

TenMomentumInt::usage = "";

EightMomentumInt::usage = "";

SixMomentumInt::usage = "";

FourMomentumInt::usage = "";

TwoMomentumInt::usage = "";

LorentzDecompose::usage = "LorentzDecompose[expr, var] decomposes all products of Scal[var,p1_MomentumQ]*...*Scal[var,pi_MomentumQ] (up to i equal to 8) and writes them like Scal[var, index1] Scal[p1, index1] ... . It can handle only at most 4 different momenta.";

TaylorExpansion::usage = "TaylorExpansion Taylor expands denominators of the form
AD[l___, den[q+k_, m_], r___],where q is a loop momentum or the sum of loop momenta
in the external momenta k up to  second order. Note that loop momenta have to be declared withDeclareLoopMomentum[expr]";

TaylorMass::usage = "TaylorMass expands denominators of the form
AD[l___, den[q_, m_], r___],where q is a loop  momentum,  in  m up to  second order, if m is NOT declared with DeclareHeavyMass[expr] as heavy mass.";                   
                    
                    
                    
TensorTwo::usage = "
TensorTwo[expr, var,var2] performs a two dimensional tensor reduction of expressions expr with up to five Lorentz Indices
assuming that the denominator of expr is an arbitrary scalar function of the variables  var1  and  var2 (Note that factorising two-loop integrals have to be tensor reduced with TensorOne).
If the numerator of  expr  depends only on  var  (var2) it performs a one-dimensional tensor reduction  in  var  (var2) using TensorOne[expr,var]  (TensorOne[expr,var2]).
Expressions like  Scal[var,q1] are treated like Scal[var,lor1]Scal[var,lor1], where  lor1 is a Lorentz index. This artificially increases the number of used Lorentz indices.Therefore it is recommended to set all quadratic scalar products to a dummy variable before performing the tensor reduction. 
.";




TensorTwoEps::usage = "TensorTwoEps[expr,var,var2]
works like TensorTwo without giving out the various intermediate
messages that indicate what the program is doing at a given time 
and setting all terms with eps^n with n>2 are set to zero."


Substitution::usage="Substitution[expr] makes substitution in the integrands of factorizing two-loop-integrals such that the propagator structure contains no overlapping loop momenta.";


TensorOne::usage = "
TensorOne[expr, var] performs the one-dimensional tensor reduction in var.
It assumes that the denominator of  expr  is an arbitrary scalar function depending on
Lorentz invariants of var. It can handle expressions $expr$ with up to 9 Lorentz Indices.
";


TensorOne6::usage = "TensorOne6[expr,var] works just like TensorOne, with one difference: it kicks out all terms proportional to var^6 or higher. This can come in handy when those terms are above your oder considered and you want to save time."



TensorOneWithoutOdd::usage = "Like TensorOne withoud dumping the odd powers of var .";


CompleteTheSquare::usage = "";

FeynParamAndMore::usage = "FeynParamAndMore[a,var]: The expression a is assumed to be of the  form  {den1, pow1, fp1}, {den2, pow2, fp2},  ... , {denN, powN, fpN}, which should mean: 1/den^pow1 1/den2^pow2*...*1/denN^powN. The variables fpI are the desired Feynman parameters. FeynParamAndMore[a,var] scrams the denominators together to form a denominator of the form (var^2-FeynDenom[1])^FeynPow[1]. It returns both FeynDenom[1] and FeynPow[1] and the summand you need to add to var in the
numerator to complete the square: AddToCompleteSquare[1]. Also it returns the prefactor FeynCoeff[1].";

FeynParamAndMore2::usage = "FeynParamAndMore2[a,var,number] works pretty much like its ancestor FeynParamAndMore. With the extra parameter, you can save the output-variables e.g. in FeynPow[number, 1] etc instead of FeynPow[1]."

dToEps::usage = "dToEps replaces d with '4 - 2 eps'.";

ScalToTimes::usage = "ScalToTimes replaces all occurences of Scal[a,b] with Times[a,b].";

gsToAlphas::usage = "gsToAlphas does what you'd guess... ."; 

RemoveSilent::usage="RemoveOddPowersOfIntVar without the messages.";

RemoveTooHigh::usage="RemoveTooHigh[expr,list,pow] removes all terms which contributions are of higher order in an expansion parameter than wanted. The second argument, lista, is a list containing all integration variables considered. The third parameter gives the highest order that contributes.";

FeynmanParametrization::usage="";

OrderDump::usage="OrderDump[expr,list,pow] is a generalization of RemoveTooHigh. To each variable in list one can give the corresponding power of that variable in the power counting.";

FourFeyn::usage="FourFeyn[list1,list2] is used to introduce three Feynman parameters (sited in list2) to pull the four denominators listed in list1 together.";
AlphasTogs::usage = "AlphaTogs does what you'd guess: it replaces \{Alpha]s with gs^2/(4 Pi).";

PartialFractionTwo::usage = " PartialFractionTwo[expr] makes partial fraction of the denominators in the two -loop case. Not approbriate for one - loop calculation, as it sets denominators with only one loop- momentum to 0. Furthermore it simplifies structures like var^2/ propagatordenominator to simpler structures 1/denominator. ";

PartialFractionOne::usage = " Partial Fraction[expr] makes partial fraction of the denominators in the one -loop case.  Furthermore it simplifies structures like var^2/ propagatordenominator to simpler structures 1/denominator. ";



DeclareHeavyMass::usage "=DeclareHeavyMass[M] declares all heavy masses. Declaration of heavy masses needed for correct usage of TaylorMass. ";

DeclareSmallMass::usage "=DeclareSmallMass[M] declares all small masses. Declaration of small masses needed for correct usage of Scaling.";

DeclareLoopMomentum::usage ="Declare any loop momentum p you wish to use with DeclareLoopMomentum[p]. Distinguishing between loop momenta and external momenta needed for correct TaylorExpansion and Scaling.";

DeclareExternalMomentum::usage =
    "Declare any external momentum p you wish to use with DeclareExternalMomentum[p]. Distinguishing between external momenta and loop  momenta needed for correct usage of TaylorExpansion and Scaling.";     

xrule::usage ="All terms x^n with n>2 are set to 0."
Scaling::usage = "All external momenta (declared with DeclareExternalMomentum) and all masses (declared with DeclareSmallMass) are multiplied with a factor x. All powers x^n with n>2 are set to zero.";

LoopMomentumQ::usage = "";
NotLoopMomentumQ::usage = "";

ExternalMomentumQ::usage = "";
NotExternalMomentumQ::usage = "";

HeavyMassQ::usage = "";
NotHeaHeavyMassQ::usage = "";

SmallMassQ::usage = "";
NotSmallMassQ::usage = "";


ClausenCl::usage =
"ClausenCl[n,x] gives the Clausen function of order n."

ClausenCl /: MakeBoxes[ClausenCl[n_,z_],TraditionalForm]:=
	RowBox[{SubscriptBox["Cl", MakeBoxes[n,TraditionalForm]],
	 "(", " ", MakeBoxes[z,TraditionalForm], ")"}]

epsrule::usage="Sets all terms with higher powers than eps^2 to 0.";

nerules::usage="Allows to express prefactors Ne[M1]Ne[M2] in terms proportional to Ne[M2], if M2 is expressed as M2=Sqrt[x1]*M1."



  

(******************************************************************************************************)

(*********************************** Start the Package ***********************************************)
Begin["`Private`"]







Off[General::spell, General::spell1]
Off[Syntax::"stresc"]



(******************************************************************************
                                                Constants
******************************************************************************)
(*Define the constants needed to shorten the result of the two-loop boxes and penguins  in the SM and THDM case*)


  (*masses and momenta*)
  (*  Print["Before usage of the functions of Integrals, you have to define which type of diagramm you want to calculate by setting values to flag (define quark type in the loop), model (defining if SM or THDM diagrams should be calculated). "]

   *)
  

  
  (*W - BOXES and Z - penguins*)
  
  
 
  
constantsWBOX:={Global`MC->0,Global`MS->0,Global`ME->0, Global`MB->0,Global`MU->0,Global`MD->0, 
			 Global`k2:>0,Global`k1:>
		     0, Global`p1:> 0, Global`p2:> 0} ;


(*SM -Modell, TOP*)
  constantsTOP:={Global`MC->0,Global`MU->0, Global`MD->0,Global`ME->0,Global`MS->0,Global`MHp->Global`MH} ;

(*SM -Modell, Charm*)
  constantsCHARM:={Global`MC->0,Global`MU->0,Global`MD->0,Global`ME->0,Global`MS->0,Global`MT->0} ;

      (*SM -model, Selbstenergie*)
	constantsselbst:={Global`MC->0,Global`MU->0,Global`MD->0,Global`ME->0} ;

Global`constants := Which[Not[FreeQ[Global`flag,Global`TOP]],constantsTOP,
			  Not[FreeQ[Global`flag,Global`CHARM]], constantsCHARM,
			  Not[FreeQ[Global`flag,Global`Z ] ],constantsWBOX,
			  Not[FreeQ[Global`flag,Global`gluon ] ],constantsTOP,
			  Not[FreeQ[Global`flag,Global`SELBST  ] ],constantsselbst]
 
Global`constantsM1M2:=
  Which[Not[FreeQ[Global`model,Global`SM]],{Global`MW->Global`M1,Global`MT->Global`M2},
	Not[FreeQ[Global`model,Global`THDM]],{Global`MH->Global`M1, Global`MT->Global`M2},
       FreeQ[Global`model,Global`SM],{Global`MW->Global`MW,Global`MT->Global`MT}]  
  
      Global`zrule:=Which[Not[FreeQ[Global`flag3,Global`gamma]],{Global`k1->Global`k1, Global`k2->Global`k2},
		    Not[FreeQ[Global`flag3,Global`Z]],{Global`k1->0,Global`k2->0,Global`MB->0}]
  
  
  
  (*****************************PREFACTORS**************************)
	
prefactorw = Global`GS^2/(16 Pi^4)* Global`EL^4/(16*Pi^4*Global`SW^4*Global`MW^2);
prefactorsm = (Global`CKM[3,3] Conjugate[Global`CKM[3,2]]*Global`GS^2*Global`EL^3)/(512*Pi^8*Global`SW^2*Global`MW^2);
prefactorgluon = (Global`CKM[3,3] Conjugate[Global`CKM[3,2]]Global`GS^2*Global`EL^2)/(512*Pi^8*Global`SW^2*Global`MW^2);
Global`prefactor:= Which[Not[FreeQ[Global`flag,Global`TOP]],prefactorsm,
		  Not[FreeQ[Global`flag,Global`CHARM]], prefactorsm,
		  Not[FreeQ[Global`flag,Global`Z ] ],prefactorsm,
			Not[FreeQ[Global`flag,Global`W ] ],prefactorw,
			Not[FreeQ[Global`flag,Global`SELBST]],prefactorsm,
             	        Not[FreeQ[Global`flag,Global`gluon]],prefactorgluon];


  
(*******************Define a new type of variable LoopMomentum which has to be declared as q1 or q2*)     

     DeclareLoopMomentum[e_] := 
    (
     DeclareMomentum[e];
     LoopMomentumQ[e]^= True;
     Unprotect[MomentumQ];
     MomentumQ[e] = True;
     Protect[MomentumQ];
     );



DeclareExternalMomentum[e_] := 
(
 DeclareMomentum[e];
     ExternalMomentumQ[e]^= True;
     Unprotect[MomentumQ];
     MomentumQ[e] = True;
     Protect[MomentumQ];
 );




DeclareHeavyMass[e_] := 
(
 DeclareMass[e];
 HeavyMassQ[e]^= True;
 Unprotect[MassQ];
 MassQ[e] = True;
 Protect[MassQ];
);

DeclareSmallMass[e_] := 
(
 DeclareMass[e];
 SmallMassQ[e]^= True;
 Unprotect[MassQ];
     MassQ[e] = True;
     Protect[MassQ];
 );




DeclareHeavyMass[p1_,p2___] :=( DeclareHeavyMass[p1] ; DeclareHeavyMass[p2]);
HeavyMassQ[c_?NumberQ p_]:= HeavyMassQ[p]
HeavyMassQ[Plus[p_,q___]] := HeavyMassQ[p] && HeavyMassQ[Plus[q]];
NotHeavyMassQ[p_]:=!HeavyMassQ[p];


DeclareSmallMass[p1_,p2___] :=( DeclareSmallMass[p1] ; DeclareSmallMass[p2]);
SmallMassQ[c_?NumberQ p_]:= SmallMassQ[p]
SmallMassQ[Plus[p_,q___]] := SmallMassQ[p] && SmallMassQ[Plus[q]];
NotSmallMassQ[p_]:=!SmallMassQ[p];


DeclareLoopMomentum[p1_,p2___] :=( DeclareLoopMomentum[p1] ; DeclareLoopMomentum[p2]);
LoopMomentumQ[c_?NumberQ p_]:= LoopMomentumQ[p]
LoopMomentumQ[Plus[p_,q___]] := LoopMomentumQ[p] && LoopMomentumQ[Plus[q]];
NotLoopMomentumQ[p_]:=!LoopMomentumQ[p];


DeclareExternalMomentum[p1_,p2___] :=( DeclareExternalMomentum[p1] ; DeclareExternalMomentum[p2]);
ExternalMomentumQ[c_?NumberQ p_]:= ExternalMomentumQ[p]
ExternalMomentumQ[Plus[p_,q___]] := ExternalMomentumQ[p] && ExternalMomentumQ[Plus[q]];
NotExternalMomentumQ[p_]:=!ExternalMomentumQ[p];



(********************* Translation of the FeynArts Output***************)






(*Translation of the FeynArts Output***************)

FeynArtsToMasterTwo[expr_] := 
Module[{res,a,b,c},
       res= expr//.{Global`FeynAmpDenominator -> AD,
       Global`FourMomentum[Global`Internal,Global`a_]->Global`q[Global`a],
      Global`FourMomentum[Global`Outgoing,Global`a_]->Global`k[Global`a],
      Global`PropagatorDenominator -> den, Global`FermionChain -> fc, 
       Global`Index[Global`Lorentz,Global`a_] -> Global`lorindex[Global`a],
Global`FeynAmpList[a___]->1;    
Global`FeynAmp[a_,b_,c___]->c,
Global`DiracSpinor[a___]->Times[],
Global`ChiralityProjector[-1]-> L,   
Global`ChiralityProjector[1]->R,
Global`DiracSlash[a_]->a,
Global`DiracMatrix[Global`lorindex[a_]]->Global`lorindex[a],
Global`MetricTensor
->Scal,
  Global`FourVector->Scal};
(*Global`Conjugate[Global`PolarizationVector][V[c_], a___,Global`lorindex[b_]]->Scal[Global`conjeps,Global`lorindex[b]]};*)
res=res//.Global`q[i_]:>ToExpression[StringJoin[ToString[q],ToString[i]]];
 
res=res//.Global`k[i_]:>ToExpression[StringJoin[ToString[k],ToString[i]]];
res=res//.Global`Index[k_, l_] :> StringJoin[ToString[k], ToString[l]];
res=res//.Global`lorindex[i_]:>
 ToExpression[StringJoin[ToString[lor],ToString[i]]];
res=res//.{lor1->Global`lor1, lor2->Global`lor2,lor3->Global`lor3,lor4->Global`lor4,lor5->Global`lor5,lor6->Global`lor6,lor7->Global`lor7,lor8->Global`lor8,lor9->Global`lor9,lor10->Global`lor10,lor11->Global`lor11};
res=res//.{k1->Global`k1, k2->Global`k2,k3->Global`k3,k4->Global`k4,k5->Global`k5,k6->Global`k6, k7->Global`k7,k8->Global`k8,k9->Global`k9,k10->Global`k10,k11->Global`k11};
 res=res//.{q1->Global`q1, q2->Global`q2,q3->Global`q3,q4->Global`q4,q5->Global`q5,q6->Global`q6, q7->Global`q7,q8->Global`q8,q9->Global`q9,q10->Global`q10,q11->Global`q11};
res=res/.Times[a___,fc[b___],c___]:>scalstructure[Times[a,c]]*fc[b];
res=res//.fc[a___, b_ + c_ + d___, e___] -> fc[a, b, e] + fc[a, c + d, e];
res=res//.Global`NonCommutative->Dirac;
res=res//.Dirac[a_]->a;
res=res//.fc[a___,Dirac[b_,c_]]->fc[a,b,c];
res=res//.Times[a___,Dirac[b_,c_],d___]->dot[a,b,c,d];
res=res//.fc[k___,dot[a___],l___]->fc[k,a,l];
res=res/.Times->hugo;
res=res//.fc[k___,hugo[a___],l___]->fc[k,a,l];
res=res//.hugo->Times;
res=res//.fc->Dirac;
 res=res//.{Dirac[l___,a_*R,r___]:>Dirac[l,a,R,r], 
        Dirac[l___,a_*L,r___]:>Dirac[l,a,L,r]};

res=res//. Dirac[a___,b_,c___]:>
    b Dirac[a,c]/;
      b=!=R&&b=!=L&&Head[b]=!=lorindex&&Head[b]=!=Plus&&MomentumQ[b]=!=True&&
        IndexQ[b]=!=True;

    res=res//.scalstructure[a___]:>a;
    res=res//DiracLinearity;
     res=res//.Dirac[a___,b_,c___]:>
    b Dirac[a,c]/;
      b=!=R&&b=!=L&&Head[b]=!=lorindex&&Head[b]=!=Plus&&MomentumQ[b]=!=True&&
        IndexQ[b]=!=True;
DeclarePolarizationVector[poleps];
res=res//.Conjugate[Global`PolarizationVector][Global`V[10], a_, b_] ...
    -> Scal[DiracAdjunction[Global`poleps],b];
res=res//.{Global`Gluon->G, Global`Colour->C};
Return[res];
];




(************************************************************************
                                          Color          
Definitions for SU(n) with n=3
************************************************************************)

  (*SUNT[a,b,c] stands for T^a_{b c}and is the color generators of SU(3)_c,SUNT[a,b] stands ...
   forT^a T^b**)
  (*SUNF[a,b,c] are the structure constants*)

  n=3; CF=4/3;





     generator:=Dispatch[
  {
 
      
      (*With Einstein's sum rule*)
      SUNT[a_,i_,k_]*SUNT[a_,k_,j_]->
       (n^2-1)/(2 n) KroneckerDelta[i,j],
      SUNT[a_,i_,j_]*SUNT[a_,k_,l_]->
       1/2 Global`Kro[i,l]*KroneckerDelta[j,k]-
       1/(6)*
       KroneckerDelta[i,j]*
       KroneckerDelta[k,l]
  }
];


colorrules:= Dispatch[
     {
       SUNF[b_,a_,c_]*SUNT[c_,b_]
       ->1/2 *I*n*SUNT[a],
       SUNF[a_,b_,c_]*SUNT[c_,b_]
       ->-1/2 *I*n*SUNT[a],
       SUNF[a_,b_,c_]*SUNT[b_,c_]
       ->+1/2 *I*n*SUNT[a],
       SUNT[k___,b_,b_,j___]->CF SUNT[k,j],
       SUNT[j___,a_,b_,a_,k___]->-1/6 SUNT[j,b,k],
       SUNT[c_,d_] SUNF[d_,b_,a_]SUNF[a_,c_,b_]->4,
       SUNT[a_,e_]SUNF[a_,d_,c_]SUNF[d_,e_,k_]SUNF[c_,k_,b]->- I 9/4 SUNT[b],
       SUNT[a_,c_,d_]SUNF[a_,d_,e_]SUNF[e_,c_,b_]->0
     }
     ];

structureconstants:=Dispatch[
    {SUNF[a_,c_,d_]*SUNF[b_,c_,d_]->n KroneckerDelta[a,b],
     SUNF[c_,b_,a_]*SUNF[d_,e_,b_]  SUNF[k_,a_,e_]->-3/2 SUNF[c,d,k],
     SUNF[a_,b_,e_]*SUNF[e_,c_,d_]+SUNF[c_,b_,e_]*SUNF[a_,e_,d_]+ 
     SUNF[d_,b_,e_]*SUNF[a_,c_,e_]->0}
];







Color[expr__] := 
  Module[{part},
	part = expr; 
         part = part//.generator;
	 part = part//.colorrules;
         part = part //.structureconstants;
         part=part//.Global`SumOver[a___]->1;
	Return[part];]



(*****************************)
(*                end color          *)
(*****************************)



(*********************************************************************
Taylor Expansion
*********************************************************************)


taylormass := AD[l___, den[q1_, M2_], r___] :> 
   AD[l, den[q1, 0], r] + AD[l, den[q1, 0], den[q1, 0], r]*
      M2^2 /; FreeQ[M2, _?HeavyMassQ] && LoopMomentumQ[q1] && 
     FreeQ[M2, 0]
TaylorMass[expr__] := Module[{part, part1, part2, part3, xx, 
    pippo}, Clear[pippo, pippo2]; part = expr; 
    part = part //. taylormass; Null; Return[part]; ]	 
     

taylorlist2 := AD[l___, den[(q1_) + (k1_), m_], r___] :> 
   AD[l, den[q1, m], r] - AD[l, den[q1, m], den[q1, m], r]*Scal[k1, k1] - 
     2*AD[l, den[q1, m], den[q1, m], r]*Scal[k1, q1] + 
     4*AD[l, den[q1, m], den[q1, m], den[q1, m], r]*Scal[k1, q1]^2 /; 
			     FreeQ[k1,_?LoopMomentumQ] && LoopMomentumQ[q1]
TaylorExpansion[expr__] := 
  Module[{part, part1, part2, part3,xx,pippo},
 part = expr; 
part=part//.taylorlist2;
Return[part];]


HeavyMassExpansion[expr__] := 
  Module[{part, part1, part2, part3,xx,pippo},
 part = expr; 
part=TaylorMass[part];

part=TaylorExpansion[part];

Return[part];]




(*zrule = {k1 -> 0, k2 -> 0};*)




(*************************************************************************************
                                           Scaling
*****************************************************************************************************************)


xrule:=Global`x^n_:>0/;n>2;

scalerule :={p_?ExternalMomentumQ:> Global`x*p, m_?SmallMassQ :> Global`x*m}




ScaleFunction[expr__] := 
     Module[{part, part1, part2, part3,xx,pippo},
	part = expr; 
part=part/.scalerule;
part=part//DiracLinearity;
part = part //. xrule;
Return[part];]
     
     
    Scaling[expr__] := 
     Module[{part, part1, part2, part3,xx,pippo},
	part = expr; 
part = Map[ScaleFunction,part];
part=part//ExpandAll;
	 part=part//.xrule;
Return[part];]
     
     


(***************************************************************************************************************)

(********************************************TensorReduction*************************************************)             

(******************************************************************************************************************)


(* Remove odd powers of integrands *)
RemoveOddPowersOfIntVar[expr_,var_] := 
    Module[{xx,ex2,pow}, 
	  ex2 = expr/.{var -> xx var};
	  ex2 = ex2 //DiracAlgebra;
	  ex2 = Expand[ex2];
	  pow = Exponent[ex2,xx];
	  Print["Removed odd powers up to ",var,"^",pow];
	  ex2 = ex2/.{xx^n_ :> 1 /; EvenQ[n]};
	  ex2 = ex2/.{xx->0};
	  Return[ex2];
      ];

RemoveTooHigh[expr_,lista_,pow_] := 
    Module[{xx,ex2}, 
      ex2=expr;
	  Do[ex2 = ex2 /. {Part[lista,i] -> Part[lista,i]*xx},{i,1,Length[lista]}];
	  ex2=ex2//DiracAlgebra;
	  ex2 = ex2 /. {xx^n_ :> 0 /; Re[n] > pow};
	  ex2 = ex2 /. {xx->1};
	  Return[ex2];
      ];

OrderDump[expr_,lista_,pow_]:=
    Module[{xx,ex2}, 
      ex2=expr;
	  Do[ex2 = ex2 /. {Part[lista,i,1] -> Part[lista,i,1]*xx^Part[lista,i,2]},{i,1,Length[lista]}];
	  ex2=ex2 //DiracAlgebra;
	  ex2 = ex2 /. {xx^n_ :> 1 /; Re[n] < pow+1};
	  ex2 = ex2 /. {xx->0};
	  Return[ex2];
      ];

RemoveSilent[expr_,var_]:=Module[{xx,ex2,pow}, 
	  ex2 = expr/.{var -> xx var};
	  ex2 = ex2 //DiracAlgebra;
	  ex2 = Expand[ex2];
	  pow = Exponent[ex2,xx];
	  ex2 = ex2/.{xx^n_ :> 1 /; EvenQ[n]};
	  ex2 = ex2/.{xx->0};
	  Return[ex2];
      ];


LorentzDecompose[expr_,r_] := 
    Module[{ex2,s=0,t=0,u=0,q=0,p=0},
	  ex2 = Expand[expr];
	  ex2 = ex2 //. {Scal[r, r] -> r^2};
	  DeclareIndex[j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,j16,j17,j18,j19,j20];
	  DeclareIndex[b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20];
	  ex2 = ex2 /.{Dirac[a___,r,b___]-> Scal[r,j1] Dirac[a,j1,b]};
	  
	  
	     (* r^8 *)
	  ex2 = ex2 /. {
	   Scal[p_?MomentumQ,r]^8 ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[p,b4] Scal[r,b5] Scal[p,b5] Scal[r,b6] Scal[p,b6]*
	   Scal[r,b7] Scal[p,b7] Scal[r,b8] Scal[p,b8],
	   
	   Scal[p_?MomentumQ,r]^7 Scal[q_?MomentumQ,r] ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[p,b4] Scal[r,b5] Scal[p,b5] Scal[r,b6] Scal[p,b6]*
	   Scal[r,b7] Scal[p,b7] Scal[r,b8] Scal[q,b8],
	   
	   Scal[p_?MomentumQ,r]^6 Scal[q_?MomentumQ]^2 ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[p,b4] Scal[r,b5] Scal[p,b5] Scal[r,b6] Scal[p,b6]*
	   Scal[r,b7] Scal[q,b7] Scal[r,b8] Scal[q,b8],
	   
	   Scal[p_?MomentumQ,r]^6 Scal[q_?MomentumQ,r] Scal[t_?MomentumQ,r] ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[p,b4] Scal[r,b5] Scal[p,b5] Scal[r,b6] Scal[p,b6]*
	   Scal[r,b7] Scal[q,b7] Scal[r,b8] Scal[t,b8],
	   
	   Scal[p_?MomentumQ,r]^5 Scal[q_?MomentumQ]^3 ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[p,b4] Scal[r,b5] Scal[p,b5] Scal[r,b6] Scal[q,b6]*
	   Scal[r,b7] Scal[q,b7] Scal[r,b8] Scal[q,b8],
	   
	   Scal[p_?MomentumQ,r]^5 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r] ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[p,b4] Scal[r,b5] Scal[p,b5] Scal[r,b6] Scal[q,b6]*
	   Scal[r,b7] Scal[q,b7] Scal[r,b8] Scal[t,b8],
	   
	   Scal[p_?MomentumQ,r]^5 Scal[q_?MomentumQ,r] Scal[t_?MomentumQ,r] Scal[s_?MomentumQ,r] ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[p,b4] Scal[r,b5] Scal[p,b5] Scal[r,b6] Scal[s,b6]*
	   Scal[r,b7] Scal[q,b7] Scal[r,b8] Scal[t,b8],
	   
	   Scal[p_?MomentumQ,r]^4 Scal[q_?MomentumQ]^4 ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[p,b4] Scal[r,b5] Scal[q,b5] Scal[r,b6] Scal[q,b6]*
	   Scal[r,b7] Scal[q,b7] Scal[r,b8] Scal[q,b8],
	   
	   Scal[p_?MomentumQ,r]^4 Scal[q_?MomentumQ,r]^3 Scal[t_?MomentumQ,r] ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[p,b4] Scal[r,b5] Scal[q,b5] Scal[r,b6] Scal[q,b6]*
	   Scal[r,b7] Scal[q,b7] Scal[r,b8] Scal[t,b8],
	   
	   Scal[p_?MomentumQ,r]^4 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r]^2 ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[p,b4] Scal[r,b5] Scal[q,b5] Scal[r,b6] Scal[q,b6]*
	   Scal[r,b7] Scal[t,b7] Scal[r,b8] Scal[t,b8],
	   	   
	   Scal[p_?MomentumQ,r]^4 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r] Scal[s_?MomentumQ,r] ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[p,b4] Scal[r,b5] Scal[q,b5] Scal[r,b6] Scal[q,b6]*
	   Scal[r,b7] Scal[s,b7] Scal[r,b8] Scal[t,b8],	 
	   
	   Scal[p_?MomentumQ,r]^3 Scal[q_?MomentumQ,r]^3 Scal[t_?MomentumQ,r]^2 ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[q,b4] Scal[r,b5] Scal[q,b5] Scal[r,b6] Scal[q,b6]*
	   Scal[r,b7] Scal[t,b7] Scal[r,b8] Scal[t,b8],
	   	   
	   Scal[p_?MomentumQ,r]^3 Scal[q_?MomentumQ,r]^3 Scal[t_?MomentumQ,r] Scal[s_?MomentumQ,r] ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[q,b4] Scal[r,b5] Scal[q,b5] Scal[r,b6] Scal[q,b6]*
	   Scal[r,b7] Scal[s,b7] Scal[r,b8] Scal[t,b8],
	   
	   Scal[p_?MomentumQ,r]^3 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r]^2 Scal[s_?MomentumQ,r] ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[p,b3]*
	   Scal[r,b4] Scal[q,b4] Scal[r,b5] Scal[q,b5] Scal[r,b6] Scal[t,b6]*
	   Scal[r,b7] Scal[t,b7] Scal[r,b8] Scal[s,b8],
	   
	   Scal[p_?MomentumQ,r]^2 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r]^2 Scal[s_?MomentumQ,r]^2 ->
	   Scal[r,b1] Scal[p,b1] Scal[r,b2] Scal[p,b2] Scal[r,b3] Scal[q,b3]*
	   Scal[r,b4] Scal[q,b4] Scal[r,b5] Scal[t,b5] Scal[r,b6] Scal[t,b6]*
	   Scal[r,b7] Scal[s,b7] Scal[r,b8] Scal[s,b8] };
	   
	  
	     (* r^7 *)
	  ex2 = ex2 /. {
	   Scal[p_?MomentumQ,r]^7 ->
	   Scal[r,b15] Scal[p,b15] Scal[r,b9] Scal[p,b9] Scal[r,b10] Scal[p,b10]*
	   Scal[r,b11] Scal[p,b11] Scal[r,b12] Scal[p,b12] Scal[r,b13] Scal[p,b13]*
	   Scal[r,b14] Scal[p,b14],
	   
	   Scal[p_?MomentumQ,r]^6 Scal[q_?MomentumQ,r] ->
	   Scal[r,b15] Scal[p,b15] Scal[r,b9] Scal[p,b9] Scal[r,b10] Scal[p,b10]*
	   Scal[r,b11] Scal[p,b11] Scal[r,b12] Scal[p,b12] Scal[r,b13] Scal[p,b13]*
	   Scal[r,b14] Scal[q,b14],
	   
	   Scal[p_?MomentumQ,r]^5 Scal[q_?MomentumQ,r]^2 ->
	   Scal[r,b15] Scal[p,b15] Scal[r,b9] Scal[p,b9] Scal[r,b10] Scal[p,b10]*
	   Scal[r,b11] Scal[p,b11] Scal[r,b12] Scal[p,b12] Scal[r,b13] Scal[q,b13]*
	   Scal[r,b14] Scal[q,b14],
	   	   
	   Scal[p_?MomentumQ,r]^5 Scal[q_?MomentumQ,r] Scal[t_?MomentumQ,r] ->
	   Scal[r,b15] Scal[p,b15] Scal[r,b9] Scal[p,b9] Scal[r,b10] Scal[p,b10]*
	   Scal[r,b11] Scal[p,b11] Scal[r,b12] Scal[p,b12] Scal[r,b13] Scal[q,b13]*
	   Scal[r,b14] Scal[t,b14],
	   
	   Scal[p_?MomentumQ,r]^4 Scal[q_?MomentumQ,r]^3 ->
	   Scal[r,b15] Scal[p,b15] Scal[r,b9] Scal[p,b9] Scal[r,b10] Scal[p,b10]*
	   Scal[r,b11] Scal[p,b11] Scal[r,b12] Scal[q,b12] Scal[r,b13] Scal[q,b13]*
	   Scal[r,b14] Scal[q,b14],
	   
	   Scal[p_?MomentumQ,r]^4 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r] ->
	   Scal[r,b15] Scal[p,b15] Scal[r,b9] Scal[p,b9] Scal[r,b10] Scal[p,b10]*
	   Scal[r,b11] Scal[p,b11] Scal[r,b12] Scal[q,b12] Scal[r,b13] Scal[q,b13]*
	   Scal[r,b14] Scal[t,b14],
	   
	   Scal[p_?MomentumQ,r]^4 Scal[q_?MomentumQ,r] Scal[t_?MomentumQ,r] Scal[s_?MomentumQ,r] ->
	   Scal[r,b15] Scal[p,b15] Scal[r,b9] Scal[p,b9] Scal[r,b10] Scal[p,b10]*
	   Scal[r,b11] Scal[p,b11] Scal[r,b12] Scal[s,b12] Scal[r,b13] Scal[q,b13]*
	   Scal[r,b14] Scal[t,b14], 
	   
	   Scal[p_?MomentumQ,r]^3 Scal[q_?MomentumQ,r]^3 Scal[t_?MomentumQ,r] ->
	   Scal[r,b15] Scal[p,b15] Scal[r,b9] Scal[p,b9] Scal[r,b10] Scal[p,b10]*
	   Scal[r,b11] Scal[q,b11] Scal[r,b12] Scal[q,b12] Scal[r,b13] Scal[q,b13]*
	   Scal[r,b14] Scal[t,b14],
	   
	   Scal[p_?MomentumQ,r]^3 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r]^2 ->
	   Scal[r,b15] Scal[p,b15] Scal[r,b9] Scal[p,b9] Scal[r,b10] Scal[p,b10]*
	   Scal[r,b11] Scal[q,b11] Scal[r,b12] Scal[q,b12] Scal[r,b13] Scal[t,b13]*
	   Scal[r,b14] Scal[t,b14],
	   
	   Scal[p_?MomentumQ,r]^3 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r] Scal[s_?MomentumQ,r] ->
	   Scal[r,b15] Scal[p,b15] Scal[r,b9] Scal[p,b9] Scal[r,b10] Scal[p,b10]*
	   Scal[r,b11] Scal[q,b11] Scal[r,b12] Scal[q,b12] Scal[r,b13] Scal[t,b13]*
	   Scal[r,b14] Scal[s,b14],
	   
	   Scal[p_?MomentumQ,r]^2 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r]^2 Scal[s_?MomentumQ,r] ->
	   Scal[r,b15] Scal[p,b15] Scal[r,b9] Scal[p,b9] Scal[r,b10] Scal[q,b10]*
	   Scal[r,b11] Scal[q,b11] Scal[r,b12] Scal[t,b12] Scal[r,b13] Scal[t,b13]*
	   Scal[r,b14] Scal[s,b14] };
	   
 
	   
	      (* r^6 *)
	  ex2 = ex2/. {
	   Scal[p_?MomentumQ,r]^6 -> 
	   Scal[r,j12] Scal[p,j12] Scal[r,j13] Scal[p,j13] Scal[r,j14] Scal[p,j14] *
	   Scal[r,j15] Scal[p,j15] Scal[r,j16] Scal[p,j16] Scal[r,j17] Scal[p,j17],
	   
	   Scal[p_?MomentumQ,r]^5 Scal[q_?MomentumQ,r] ->
	   Scal[r,j12] Scal[p,j12] Scal[r,j13] Scal[p,j13] Scal[r,j14] Scal[p,j14] *
	   Scal[r,j15] Scal[p,j15] Scal[r,j16] Scal[p,j16] Scal[r,j17] Scal[q,j17],
	   
	   Scal[p_?MomentumQ,r]^4 Scal[q_?MomentumQ,r]^2 ->
	   Scal[r,j12] Scal[p,j12] Scal[r,j13] Scal[p,j13] Scal[r,j14] Scal[p,j14] *
	   Scal[r,j15] Scal[p,j15] Scal[r,j16] Scal[q,j16] Scal[r,j17] Scal[q,j17],
	   
	   Scal[p_?MomentumQ,r]^4 Scal[q_?MomentumQ,r] Scal[t_?MomentumQ,r] ->
	   Scal[r,j12] Scal[p,j12] Scal[r,j13] Scal[p,j13] Scal[r,j14] Scal[p,j14] *
	   Scal[r,j15] Scal[p,j15] Scal[r,j16] Scal[q,j16] Scal[r,j17] Scal[t,j17],
	   
	   
	   Scal[p_?MomentumQ,r]^3 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r] ->
	   Scal[r,j12] Scal[p,j12] Scal[r,j13] Scal[p,j13] Scal[r,j14] Scal[p,j14] *
	   Scal[r,j15] Scal[q,j15] Scal[r,j16] Scal[q,j16] Scal[r,j17] Scal[t,j17],
	   
	   Scal[p_?MomentumQ,r]^3 Scal[q_?MomentumQ,r]^3 ->
	   Scal[r,j12] Scal[p,j12] Scal[r,j13] Scal[p,j13] Scal[r,j14] Scal[p,j14] *
	   Scal[r,j15] Scal[q,j15] Scal[r,j16] Scal[q,j16] Scal[r,j17] Scal[q,j17],
	   
	   Scal[p_?MomentumQ,r]^3 Scal[q_?MomentumQ,r] Scal[t_?MomentumQ,r]*
	   Scal[s_?MomentumQ,r] ->
	   Scal[r,j12] Scal[p,j12] Scal[r,j13] Scal[p,j13] Scal[r,j14] Scal[p,j14] *
	   Scal[r,j15] Scal[q,j15] Scal[r,j16] Scal[t,j16] Scal[r,j17] Scal[s,j17],
	   
	   Scal[p_?MomentumQ,r]^2 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r]^2 ->
	   Scal[r,j12] Scal[p,j12] Scal[r,j13] Scal[p,j13] Scal[r,j14] Scal[q,j14] *
	   Scal[r,j15] Scal[q,j15] Scal[r,j16] Scal[t,j16] Scal[r,j17] Scal[t,j17],
	   
	   Scal[p_?MomentumQ,r]^2 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r]*
	   Scal[s_?MomentumQ,r] ->
	   Scal[r,j12] Scal[p,j12] Scal[r,j13] Scal[p,j13] Scal[r,j14] Scal[q,j14] *
	   Scal[r,j15] Scal[q,j15] Scal[r,j16] Scal[t,j16] Scal[r,j17] Scal[s,j17] };
	   
	      (* r^5 *)
	  ex2 = ex2/. {
	   Scal[p_?MomentumQ,r]^5 -> 
	   Scal[r,b16] Scal[p,b16] Scal[r,b17] Scal[p,b17] Scal[r,b18] Scal[p,b18] *
	   Scal[r,b19] Scal[p,b19] Scal[r,b20] Scal[p,b20],
	   
	   Scal[p_?MomentumQ,r]^4 Scal[q_?MomentumQ,r] -> 
	   Scal[r,b16] Scal[p,b16] Scal[r,b17] Scal[p,b17] Scal[r,b18] Scal[p,b18] *
	   Scal[r,b19] Scal[p,b19] Scal[r,b20] Scal[q,b20],
	   
	   Scal[p_?MomentumQ,r]^3 Scal[q_?MomentumQ,r]^2 -> 
	   Scal[r,b16] Scal[p,b16] Scal[r,b17] Scal[p,b17] Scal[r,b18] Scal[p,b18] *
	   Scal[r,b19] Scal[q,b19] Scal[r,b20] Scal[q,b20],
	   
	   Scal[p_?MomentumQ,r]^3 Scal[q_?MomentumQ,r] Scal[t_?MomentumQ,r] -> 
	   Scal[r,b16] Scal[p,b16] Scal[r,b17] Scal[p,b17] Scal[r,b18] Scal[p,b18] *
	   Scal[r,b19] Scal[q,b19] Scal[r,b20] Scal[t,b20],
	   
	   Scal[p_?MomentumQ,r]^2 Scal[q_?MomentumQ,r]^2 Scal[t_?MomentumQ,r] -> 
	   Scal[r,b16] Scal[p,b16] Scal[r,b17] Scal[p,b17] Scal[r,b18] Scal[q,b18] *
	   Scal[r,b19] Scal[q,b19] Scal[r,b20] Scal[t,b20],
	   
	   Scal[p_?MomentumQ,r]^2 Scal[q_?MomentumQ,r] Scal[t_?MomentumQ,r] Scal[s_?MomentumQ,r] -> 
	   Scal[r,b16] Scal[p,b16] Scal[r,b17] Scal[p,b17] Scal[r,b18] Scal[q,b18] *
	   Scal[r,b19] Scal[t,b19] Scal[r,b20] Scal[s,b20] };
	   	           
	  
          (* r^4 *)
	  ex2 = ex2/. {
	      Scal[p_?MomentumQ, r]^4  ->  
		  Scal[r, j2] Scal[p, j2]  Scal[r, j3] Scal[p, j3]*
	          Scal[r, j4] Scal[p, j4] Scal[r, j5] Scal[p, j5],   
         
	      Scal[p_?MomentumQ, r]^2 Scal[q_?MomentumQ, r]^2  ->   
		  Scal[r, j2] Scal[p, j2]  Scal[r, j3] Scal[p, j3]*
	          Scal[r, j4] Scal[q, j4] Scal[r, j5] Scal[q, j5], 

	      Scal[p_?MomentumQ, r]^2 Scal[a_?MomentumQ, r]*  
	      Scal[b_?MomentumQ, r]  ->  
		  Scal[r, j2] Scal[p, j2]  Scal[r, j3] Scal[p, j3]  Scal[a, j4]* 
	          Scal[r, j4] Scal[b, j5] Scal[r, j5], 

	      Scal[a_?MomentumQ, r]  Scal[b_?MomentumQ, r]* 
	      Scal[c_?MomentumQ, r]  Scal[d_?MomentumQ, r]  ->    
		  Scal[a, j2] Scal[r, j2]  Scal[b, j3] Scal[r, j3]  Scal[c, j4]*
                  Scal[r, j4] Scal[d, j5] Scal[r, j5]                          }; 


          (* r^3  (together with Scal[r, j5] Dirac[a, j5, b]) *) 
	  ex2 = ex2 /. {
	      Scal[p_?MomentumQ, r]^3  ->    
		  Scal[r, j6] Scal[p, j6]  Scal[r, j7] Scal[p, j7] Scal[r, j8] Scal[p, j8], 
	      Scal[p_?MomentumQ, r]^2 Scal[q_?MomentumQ, r]  ->   
		  Scal[r, j6] Scal[p, j6]  Scal[r, j7] Scal[p, j7] Scal[r, j8] Scal[q, j8],
		  Scal[a_?MomentumQ, r]  Scal[b_?MomentumQ, r] Scal[c_?MomentumQ, r]   ->    
		  Scal[a, j6] Scal[r, j6]  Scal[b, j7] Scal[r, j7]  Scal[c, j8] Scal[r, j8]  }; 

          (* r^2 *) 
	  ex2 = ex2 /. {
	      Scal[a_?MomentumQ, r]^2 ->  
		  Scal[r, j9] Scal[a, j9] Scal[a, j10] Scal[r, j10],		  
		  Scal[a_?MomentumQ, r] Scal[b_?MomentumQ, r]  ->  
		  Scal[r, j9] Scal[a, j9] Scal[b, j10] Scal[r, j10] }; 


          (* r (together with Scal[r,j5]Dirac[a,j5,b]) *)
	  ex2 = ex2/.{Scal[p_?MomentumQ,r]-> Scal[p,j11] Scal[r,j11]};
      ex2 = ex2 //. {r^n_?EvenQ :> Scal[r,r]^(n/2)};
	  Return[ex2];

];

FourScal[a_,b_,c_,e_] := Scal[a,b] Scal[c,e] + Scal[a,c] Scal[b,e] + Scal[a,e] Scal[b,c];

SixScal[al_,be_,ga_,de_,rh_,si_]:= Scal[al,be]*(Scal[de,ga] Scal[rh,si] + Scal[de,rh] Scal[ga,si] + Scal[de,si] Scal[ga,rh]) +
                       Scal[al,de]*(Scal[be,ga] Scal[rh,si] + Scal[be,rh] Scal[ga,si] + Scal[be,si] Scal[ga,rh]) +
                       Scal[al,ga]*(Scal[be,de] Scal[rh,si] + Scal[be,rh] Scal[de,si] + Scal[be,si] Scal[rh,de]) +
                       Scal[al,rh]*(Scal[be,ga] Scal[de,si] + Scal[be,de] Scal[ga,si] + Scal[be,si] Scal[ga,de]) +
                       Scal[al,si]*(Scal[be,ga] Scal[de,rh] + Scal[be,de] Scal[ga,rh] + Scal[be,rh] Scal[ga,de]);
 
EightScal[a_,b_,c_,e_,f_,g_,h_,k_]:= Expand[Scal[a,b] SixScal[c,e,f,g,h,k] + Scal[a,c] SixScal[b,e,f,g,h,k] +
                       Scal[a,e] SixScal[c,b,f,g,h,k] + Scal[a,f] SixScal[b,e,c,g,h,k] +
                       Scal[a,g] SixScal[c,b,f,e,h,k] + Scal[a,h] SixScal[b,e,c,g,f,k] +
                       Scal[a,k] SixScal[b,c,e,f,g,h]];
TenScal[a_,b_,c_,e_,f_,g_,h_,k_,l_,m_] := Expand[Scal[a,b] EightScal[c,e,f,g,h,k,l,m] + Scal[a,c] EightScal[b,e,f,g,h,k,l,m] +
                       Scal[a,e] EightScal[c,b,f,g,h,k,l,m] + Scal[a,f] EightScal[b,e,c,g,h,k,l,m] +
                       Scal[a,g] EightScal[c,e,f,b,h,k,l,m] + Scal[a,h] EightScal[b,e,f,g,c,k,l,m] +
                       Scal[a,k] EightScal[c,e,f,g,h,b,l,m] + Scal[a,l] EightScal[b,e,f,g,h,k,c,m] +
                       Scal[a,m] EightScal[b,c,e,f,g,h,k,l]];                       





	TwoScal[a_,b_,c_,e_] :=  Scal[a,c] Scal[b,e] + Scal[a,e] Scal[b,c];


	TwoMomentumTwoInt[r_,k_]:= {Scal[a_?IndexQ, r] Scal[b_?IndexQ, k] -> 
        (1/4 + eps/8 + eps^2/16)  Scal[r, k] Scal[a, b]};

	FourMomentumTwoInt[r_,k_] := {Scal[a_, r] Scal[b_, r] Scal[c_, r] Scal[e_, k] -> 
       (1/24+5*eps/144+ 19 *eps^2/864)   (Scal[r, r]* Scal[r,k]*FourScal[a, b, c, e])};



	FourMomentumTwoIntTwo[r_,k_] := {Scal[a_, r] Scal[b_, r] Scal[c_, k] Scal[e_, k] -> 
					 (1/72+eps/48+55 *eps^2/2592)
(((1+(4-2*eps))Scal[r, r]* Scal[k,k]- 2*Scal[r,k]^2)* Scal[a,b]*Scal[c,e]+
				    ((4-2*eps)*Scal[r,k]^2-Scal[r,r]*Scal[k,k]
                                     )*TwoScal[a, b, c, e])};









	(*1/d is replaced by   (1/4 +eps/8+ eps^2/16) *)
TwoMomentumInt[r_] := {Scal[a_?IndexQ, r] Scal[b_?IndexQ, r] ->
 (1/4 +eps/8+ eps^2/16) Scal[r,r] Scal[a, b] };


FourMomentumInt[r_]:=
	{Scal[a_, r] Scal[b_, r] Scal[c_,r] Scal[e_,r] ->  
        (1/24+5*eps/144+(19*eps^2)/864) Scal[r,r]^2 FourScal[a,b,c,e] };

SixMomentumInt[r_] := {Scal[al_, r] Scal[be_, r] Scal[de_, r] *
                       Scal[ga_, r] Scal[r,rh_] Scal[r, si_] -> 
                       (1/192+ (13*eps/2304)+(115*eps^2)/27648) Scal[r, r]^3 * SixScal[al,be,de,ga,rh,si] };
                       
EightMomentumInt[r_]:={Scal[a_?IndexQ, r] Scal[b_?IndexQ, r] Scal[c_?IndexQ, r] Scal[e_?IndexQ, r] *
                       Scal[f_?IndexQ, r] Scal[r, g_?IndexQ] Scal[r, h_?IndexQ] Scal[ke_?IndexQ, r] ->
                      (1/1920+77*eps/115200+3799*eps^2/6912000) Scal[r,r]^4 * EightScal[a,b,c,e,f,g,h,ke] };
                      
           
TenMomentumInt[r_] := {Scal[a_?IndexQ, r] Scal[b_?IndexQ, r] Scal[c_?IndexQ, r] Scal[e_?IndexQ, r] *
                       Scal[f_?IndexQ, r] Scal[r, g_?IndexQ] Scal[r, h_?IndexQ] Scal[k_?IndexQ, r] *
                       Scal[l_?IndexQ, r] Scal[m_?IndexQ,r] -> 
                       1/(d (2 + d) (4 + d) (6 + d) (8 + d) ) Scal[r,r]^5 * TenScal[a,b,c,e,f,g,h,k,l,m]};
                      


epsrule := eps^x_ :> 0 /; x > 2




 TensorOne[expr_,var_]:=
    Module[{res,i,res2,xx,ex2,pow,meinlist={},endpow,toohigh={}},
      Clear[ema,j,a,m]; 

       ex2=DenruleOne[expr,var];     
      ex2=ex2/.{var->xx var};
      ex2=ex2//DiracAlgebra;
      ex2=Expand[ex2];
      pow=Exponent[ex2,xx];
      ex2=ex2/.{xx^n_:>1/;EvenQ[n]};
      ex2=ex2/.{xx->0};
      res=ex2;
      res=LorentzDecompose[res,var];
      res=res//Expand;
      res2=res//. Scal[var,_?IndexQ]->ema;
      res2=res2//Expand;
      meinlist=Cases[res2,ema^_,{0,Infinity}]//Union;
      For[i=1,i<Length[meinlist]+1,i++,pott[i]=meinlist[[i,2]]];
      endpow=Max[Table[pott[i],{i,1,Length[meinlist]}],2];
      toohigh=Cases[res,Scal[var,_?IndexQ]^Condition[au_,Re[au]>1],{0,Infinity}]//Union;
      If[Length[toohigh]>0||endpow>8, Print["Sorry, but the program ...
                          can handle index structures up to ",var,"^9 ...
                          only! TensorOne aborted... "],
      If[endpow>7, res=res//.EightMomentumInt[var];];
      If[endpow>5,res=res//.SixMomentumInt[var];
      ];
      If[endpow>3, res=res//.FourMomentumInt[var];];

      res=res//.TwoMomentumInt[var];
      res=Expand[res];
      res=
        ContractIndex[
          res,{j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,j16,j17,j18,
            j19,j20}];
      res=
        ContractIndex[
          res,{b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,
            b19,b20}];
      res=res//DiracAlgebra;
      res=res/.dToEps;
     

      res=Expand[res];
      res=res//.epsrule;
     Return[DenruleOneBack[res,var]];
]];



 TensorOneWithoutOdd[expr_,var_]:=
    Module[{res,res2,xx,yy,ex2,pow,pow2,meinlist={},meinlist2={},endpow,endpow2,toohigh={},toohigh2={}},
      Clear[ema,a,m];
      Print["Works safely up to ",var,
        "^9 only! Also, only up to four different momenta are allowed!"]; 
	   res=DenruleOne[expr,var];   
      res=res//DiracAlgebra//Expand;
      res=LorentzDecompose[res,var];
      res=res//Expand;
      Print["Done with LorentzDecompose."];
      res2=res//. Scal[var,_?IndexQ]->ema;
      res2=res2//Expand;
      meinlist=Cases[res2,ema^_,{0,Infinity}]//Union;
      For[i=1,i<Length[meinlist]+1,i++,pott[i]=meinlist[[i,2]]];
      endpow=Max[Table[pott[i],{i,1,Length[meinlist]}],2];
      toohigh=Cases[res,Scal[var,_?IndexQ]^Condition[au_,Re[au]>1],{0,Infinity}]//Union;
      If[Length[toohigh]>0, Print["Sorry, but the program can handle index structures up to ",var,"^8 only! TensorOneaborted... "],
      Print[endpow];
      If[endpow>7, res=res//.EightMomentumInt[var];      
       Print["Done with EightMomentum."]];
      If[endpow>5,res=res//.SixMomentumInt[var];
       Print["Done with SixMomentum."]];
      If[endpow>3, res=res//.FourMomentumInt[var];
      Print["Done with FourMomentum."]];
      res=res//.TwoMomentumInt[var];
      Print["Done with TwoMomentum."];
      res=Expand[res];
      res=
        ContractIndex[
          res,{j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,j16,j17,j18,
            j19,j20}];
      res=
        ContractIndex[
          res,{b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,
            b19,b20}];
      res=res//DiracAlgebra;
      res=res/.dToEps//Factor;
	 Return[DenruleOneBack[res,var]];
      ]];   



TensorTwoEps[expr_,var_,var2_]:=
  Module[{res,res2,i,ex2,ex1,pow,pow1,pow2,xx,yy,meinlist={},meinlist2={},endpow,endpow2,toohigh={},toohigh2={}},
	 Clear[ema,a,j,m,var11,var12,var22];
	 Clear[ema2];
	 Clear[xx];
	 Clear[yy];
	 ex1=Denrule[expr,var,var2];
	 ex1=Dummyrule[ex1,var,var2];
	 ex2 = ex1 /. {var -> xx var, var2 -> yy var2};
	 
	 ex2=DiracLinearity[ex2];
	 
	 pow1 = Exponent[ex2, xx];
	 pow2= Exponent[ex2, yy];	 
	 If[pow1<1,
	    ex1 = DenruleBack[ex1,var, var2];
	    ex1=TensorOne[ex1,var2];
	    ex1=DummyruleBack[ex1,var,var2];
	    Return[ex1]];
	 If[pow2<1,
	    ex1 = DenruleBack[ex1, var, var2];
	    ex1=TensorOne[ex1,var];
	    ex1=DummyruleBack[ex1,var,var2];
	    Return[ex1]];
	 ex2=ex2//.yy->xx;
	 ex2=Expand[ex2];
	 pow = Exponent[ex2, xx];
	 If[pow>5,Print["Exponent in ",var,"=",pow1," ,exponent in ", var2,"=",pow2,"TensorTwo can handle terms only up to powers of 5 in the LoopMomenta. Tensor Reduction aborted."];Return[]];
	 ex2 = ex2 /. {xx^n_ :> 1 /; EvenQ[n]};
	 ex2 = ex2 /. {xx -> 0};
	 res = ex2;
	 res = LorentzDecompose[res, var];
	 DeclareIndex[u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20];
         DeclareIndex[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20];
         res = res /. {j1->u1,j2->u2,j3->u3,j4->u4,j5->u5,j6->u6,j7->u7,
	 	       j8->u8,j9->u9,j10->u10,j11->u11,j12->u12,j13->u13,
		       j14->u14,j15->u15,j16->u16,j17->u17,j18->u18,j19->u19,
		       j20->u20} /.
  {b1->v1,b2->v2,b3->v3,b4->v4,b5->v5,b6->v6,b7->v7,
   b8->v8,b9->v9,b10->v10,b11->v11,b12->v12,b13->v13,
		       b14->v14,b15->v15,b16->v16,b17->v17,b18->v18,b19->v19,
		       b20->v20};

	 res=LorentzDecompose[res,var2];
	 res = res // Expand;
	 res2 = res //. {Scal[var, _?IndexQ] -> ema, Scal[var2, _?IndexQ] -> ema2};
	 
	 res2 = res2 // Expand;
	 
	 meinlist=Cases[res2,ema^_,{0,Infinity}]//Union;
	 meinlist2 = Cases[res2, ema2^_, {0, Infinity}] // Union;
	 
	 
	 
	 For[i=1,i<Length[meinlist]+1,i++,pott[i]=meinlist[[i,2]]];
	 For[i=1,i < Length[meinlist2] + 1, i++,pott2[i] = meinlist2[[i, 2]]];
	 endpow=Max[Table[pott[i], {i, 1, Length[meinlist]}],1];
	 endpow2=Max[Table[pott2[i], {i, 1, Length[meinlist2]}],1];
	 
	 toohigh=Cases[res,Scal[var,_?IndexQ]^Condition[au_,Re[au]>1],{0,Infinity}]//Union;
	 toohigh2 = Cases[res,Scal[var2, _?IndexQ]^Condition[au_, Re[au] > 1], {0, Infinity}]//Union;
	 If[Length[toohigh]>0||Length[toohigh2]>0, 
	    Print["Sorry, but the program can handle index structures up to ",var,"^8 only! TensorOneaborted... "],
	    If[endpow>2, res=res//. FourMomentumTwoInt[var,var2];    
	       ];
	    If[endpow2>2, res=res//. FourMomentumTwoInt[var2,var];    
	       ];
	    If[endpow>1, res=res//. FourMomentumTwoIntTwo[var,var2]];
	    If[endpow >0, res = res //. TwoMomentumTwoInt[var, var2]];
	    If[endpow>7, res=res//.EightMomentumInt[var]];
	    If[endpow>5,res=res//.SixMomentumInt[var]];
	    If[endpow>3, res=res//.FourMomentumInt[var]];
	    
	     If[endpow>1, res=res//.TwoMomentumInt[var]];
		
		
		If[endpow2>7, res=res//.EightMomentumInt[var2]];
		If[endpow2>5,res=res//.SixMomentumInt[var2]];
		If[endpow2>3, res=res//.FourMomentumInt[var2]];
		
		If[endpow2>1, res=res//.TwoMomentumInt[var2]];
		   
		   res=Expand[res];
		   
		   res=
		   ContractIndex[
				 res,{j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,j16,j17,j18,
				      j19,j20}];
		   
		   res=ContractIndex[
				     res,{b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,
					  b19,b20}];
		   
		   res= ContractIndex[
				      res,{u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,
					   u19,u20}];
		   res=ContractIndex[
				     res,{v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,
			   v19,v20}];
		   
		   res=res//DiracAlgebra;
		   res=res/.dToEps//Expand;
		   res=res//.epsrule;
		   res=res//Factor;

		   res=DenruleBack[res,var,var2];
		   
		   res=DummyruleBack[res,var,var2];
		   
		   Return[res];
		   ]];   

TensorTwo[expr_,var_,var2_]:=
  Module[{res,res2,i,ex2,ex1,pow,pow1,pow2,xx,yy,meinlist={},meinlist2={},endpow,endpow2,toohigh={},toohigh2={}},
	 Clear[ema,a,j,m,var11,var12,var22];
	 Clear[ema2];
	 Clear[xx];
	 Clear[yy];
	 ex1=Denrule[expr,var,var2];
	 ex1=Dummyrule[ex1,var,var2];
	 ex2 = ex1 /. {var -> xx var, var2 -> yy var2};
	 
	 ex2=DiracLinearity[ex2];
	 
	 pow1 = Exponent[ex2, xx];
	 pow2= Exponent[ex2, yy];	 
	 If[pow1<1,
	    ex1 = DenruleBack[ex1,var, var2];
	    ex1=TensorOne[ex1,var2];
	    ex1=DummyruleBack[ex1,var,var2];
	    Return[ex1]];
	 If[pow2<1,
	    ex1 = DenruleBack[ex1, var, var2];
	    ex1=TensorOne[ex1,var];
	    ex1=DummyruleBack[ex1,var,var2];
	    Return[ex1]];
	 ex2=ex2//.yy->xx;
	 ex2=Expand[ex2];
	 pow = Exponent[ex2, xx];
	 If[pow>5,Print["Exponent in ",var,"=",pow1," ,exponent in ", var2,"=",pow2,"TensorTwo can handle terms only up to powers of 5 in the LoopMomenta. Tensor Reduction aborted."];Return[]];
	 ex2 = ex2 /. {xx^n_ :> 1 /; EvenQ[n]};
	 ex2 = ex2 /. {xx -> 0};
	 res = ex2;
	 res = LorentzDecompose[res, var];
	 DeclareIndex[u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20];
         DeclareIndex[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20];
         res = res /. {j1->u1,j2->u2,j3->u3,j4->u4,j5->u5,j6->u6,j7->u7,
	 	       j8->u8,j9->u9,j10->u10,j11->u11,j12->u12,j13->u13,
		       j14->u14,j15->u15,j16->u16,j17->u17,j18->u18,j19->u19,
		       j20->u20} /.
  {b1->v1,b2->v2,b3->v3,b4->v4,b5->v5,b6->v6,b7->v7,
   b8->v8,b9->v9,b10->v10,b11->v11,b12->v12,b13->v13,
		       b14->v14,b15->v15,b16->v16,b17->v17,b18->v18,b19->v19,
		       b20->v20};

	 res=LorentzDecompose[res,var2];
	 res = res // Expand;
	 res2 = res //. {Scal[var, _?IndexQ] -> ema, Scal[var2, _?IndexQ] -> ema2};
	 
	 res2 = res2 // Expand;
	 
	 meinlist=Cases[res2,ema^_,{0,Infinity}]//Union;
	 meinlist2 = Cases[res2, ema2^_, {0, Infinity}] // Union;
	 
	 
	 
	 For[i=1,i<Length[meinlist]+1,i++,pott[i]=meinlist[[i,2]]];
	 For[i=1,i < Length[meinlist2] + 1, i++,pott2[i] = meinlist2[[i, 2]]];
	 endpow=Max[Table[pott[i], {i, 1, Length[meinlist]}],1];
	 endpow2=Max[Table[pott2[i], {i, 1, Length[meinlist2]}],1];
	 
	 toohigh=Cases[res,Scal[var,_?IndexQ]^Condition[au_,Re[au]>1],{0,Infinity}]//Union;
	 toohigh2 = Cases[res,Scal[var2, _?IndexQ]^Condition[au_, Re[au] > 1], {0, Infinity}]//Union;
	 If[Length[toohigh]>0||Length[toohigh2]>0, 
	    Print["Sorry, but the program can handle index structures up to ",var,"^8 only! TensorOneaborted... "],
	    If[endpow>2, res=res//. FourMomentumTwoInt[var,var2];    
	       ];
	    If[endpow2>2, res=res//. FourMomentumTwoInt[var2,var];    
	       ];
	    If[endpow>1, res=res//. FourMomentumTwoIntTwo[var,var2]];
	    If[endpow >0, res = res //. TwoMomentumTwoInt[var, var2]];
	    If[endpow>7, res=res//.EightMomentumInt[var]];
	    If[endpow>5,res=res//.SixMomentumInt[var]];
	    If[endpow>3, res=res//.FourMomentumInt[var]];
	    
	     If[endpow>1, res=res//.TwoMomentumInt[var]];
		
		
		If[endpow2>7, res=res//.EightMomentumInt[var2]];
		If[endpow2>5,res=res//.SixMomentumInt[var2]];
		If[endpow2>3, res=res//.FourMomentumInt[var2]];
		
		If[endpow2>1, res=res//.TwoMomentumInt[var2]];
		   
		   res=Expand[res];
		   
		   res=
		   ContractIndex[
				 res,{j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,j16,j17,j18,
				      j19,j20}];
		   
		   res=ContractIndex[
				     res,{b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,
					  b19,b20}];
		   
		   res= ContractIndex[
				      res,{u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,
					   u19,u20}];
		   res=ContractIndex[
				     res,{v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,
			   v19,v20}];
		   
		   res=res//DiracAlgebra;
		   res=res/.dToEps//Expand;
		 
		   res=res//Factor;

		   res=DenruleBack[res,var,var2];
		   
		   res=DummyruleBack[res,var,var2];
		   
		   Return[res];
		   ]];  

	    
	    
  (* To complete a square: *)
	    (* a x^2 + b x + c == a(x + 1/2 b/a)^2 + c - 1/4 b^2/a *)
	    (* CompleteTheSquare returns what you must ADD to the variable *)
(* Works for p*q as well as for Scal[p,q] *)
CompleteTheSquare[expr_,var_] := 
    Module[{aa,bb,cc,res, tmp},
	  tmp = expr//.{Scal -> Times};
	  aa  =  Coefficient[tmp, var^2]; 
	  bb  =  Coefficient[tmp, var]; 
	  cc = tmp /. {var ->  0}; 			       
	  res = Factor[-1/2 bb/aa];
	  res = res//.{ p_?MomentumQ^2 -> Scal[p,p], p_?MomentumQ q_?MomentumQ -> Scal[p,q]};
	  Return[res];
];



    
    
    
    
    
    
    
   


DenruleOne[expr_, var_] := Module[{part,m,j},
			       part = expr /.{den[var, m_] -> den[Global`q11, m],
					       den[var + j_, m_] -> den[Global`q11 + j, m]};
    Return[part];]

DenruleOneBack[expr_, var_] := 
  Module[{part,j,m,q11}, 
       part = expr /. den[Global`q11, m_] -> den[var, m]; 
       part=part/.den[Global`q11+j_,m_]->den[var+j,m];
       Return[part]; ]

Denrule[expr_, var_, var2_] :=
Module[{part,a,m}, 
       part = expr /. 
    {den[var, m_] -> den[Global`q11, m], 
     den[var2, m_] -> den[Global`q22, m], 
     den[var + var2, m_] -> den[Global`q11 + Global`q22, m], 
     den[var - var2, m_] -> den[Global`q11 - Global`q22, m],
     den[var+a_,m_]->den[Global`q11+a,m],
     den[var2+a_,m_]->den[Global`q22+a,m]};
       Return[part]; ]
  

DenruleBack[expr_, var_, var2_] := 
  Module[{part,a,m}, 
	 part=expr/.{den[Global`q11, m_] -> den[var, m], 
		     den[Global`q22, m_] -> den[var2, m], 
		     den[Global`q11 + Global`q22, m_] -> den[var + var2, m], 
		     den[Global`q11 - Global`q22, m_] -> den[var - var2, m],
		     den[Global`q11+a_,m_]->den[var+a,m],
		     den[Global`q22+a_,m_]->den[var2+a,m]};
Return[part]; ]






Dummyrule[expr_, var_, var2_] :=
  Module[{part}, 
       part = expr /. {Scal[var, var] -> var11, 
     Scal[var,var2] -> var12, 
     Scal[var2,var2]->var22};
       Return[part]; ]
  

  DummyruleBack[expr_, var_, var2_] := 
  Module[{part}, 
	 part=expr/.{var11 -> Scal[var,var],
		     var12 -> Scal[var,var2],
		    var22->Scal[var2,var2]};
Return[part]; ]
  





(***************************End of Tensor Reduction*)


 (********************** Feynman Parametrization*********************)


FeynmanParametrization[a___] :=     
    Module[{arg, PowSum, GammaProd, denom, fpPrefac, i, res},   

	  (*  a of  the  form  {den1, pow1, fp1}, {den2, pow2, fp2},  ... , {denN, powN, fpN}  *)
	  (*  Meaning: 1/den^pow1 1/den2^pow2*...*1/denN^powN, fpI are  the desired Feynman Variables *)
	  (*  PowSum :   Sum  of  the  powers; need  that  for  Gamma[ pow1 +   ... ]  *)
	  (*  GammaProd :   the  denominator  of  the  prefactor  *)
	  (*  denom  collects  the  final  denominator  *)
	  (*  fpPrefac  collects  the  prefactors  involving  Feynamn Parameters  *)

	  arg   =   {a};  
	  For[PowSum = 0; GammaProd = 1; fpPrefac = 1; denom = 0;    
	      i = 1, i <= Length[arg],  i++,   
	      PowSum    += arg[[i, 2]];   
	      GammaProd *= Gamma[arg[[i, 2]]];   
	      fpPrefac  *= arg[[i, 3]]^(arg[[i, 2]] - 1); 
	      denom     += arg[[i, 1]] * arg[[i, 3]];
	     ]; 

	  (*   res = Gamma[PowSum]/GammaProd fpPrefac denom^(-PowSum); *)   
	      
	  res = FeynCoeff[0] FeynDenom[0]^(-FeynPow[0]);

	  FeynCoeff[1] = Gamma[PowSum]/GammaProd * fpPrefac;   
	  FeynDenom[1] = denom;
	  FeynPow[1]   = PowSum;  
       	  Print[Unevaluated[FeynCoeff[1]]," = ", FeynCoeff[1]];
	  Print[Unevaluated[FeynDenom[1]]," = ", FeynDenom[1]];
	  Print[Unevaluated[FeynPow[1]],"   = ", FeynPow[1]];
	  Return[res];

	 ];


    

FeynParamAndMore[{a___}, pr_] :=     
    Module[{arg, PowSum, GammaProd, denom, fpPrefac, i, res, Denom},   

	  (*  a of  the  form  {den1, pow1, fp1}, {den2, pow2, fp2},  ... , {denN, powN, fpN}  *)
	  (*  Meaning: 1/den^pow1 1/den2^pow2*...*1/denN^powN, fpI are  the desired Feynman Variables *)
	  (*  PowSum :   Sum  of  the  powers; need  that  for  Gamma[ pow1 +   ... ]  *)
	  (*  GammaProd :   the  denominator  of  the  prefactor  *)
	  (*  denom  collects  the  final  denominator  *)
	  (*  fpPrefac  collects  the  prefactors  involving  Feynamn Parameters  *)

	  arg   =   {a};  
	  For[PowSum = 0; GammaProd = 1; fpPrefac = 1; denom = 0;    
	      i = 1, i <= Length[arg],  i++,   
	      PowSum    += arg[[i, 2]];   
	      GammaProd *= Gamma[arg[[i, 2]]];   
	      fpPrefac  *= arg[[i, 3]]^(arg[[i, 2]] - 1); 
	      denom     += arg[[i, 1]] * arg[[i, 3]];
	     ]; 

	  (*   res = Gamma[PowSum]/GammaProd fpPrefac denom^(-PowSum); *)   
	      
	  FeynCoeff[1] = Gamma[PowSum]/GammaProd * fpPrefac;   
	  FeynDenom[1] = denom;
	  FeynPow[1]   = PowSum;  
       	  Print[Unevaluated[FeynCoeff[1]]," = ", FeynCoeff[1]];
	  Print[Unevaluated[FeynDenom[1]]," = ", FeynDenom[1]];
	  Print[Unevaluated[FeynPow[1]],"   = ", FeynPow[1]];

	   (* And Complete the Square *)
	   AddToCompleteSquare[1] = CompleteTheSquare[FeynDenom[1], pr];
	   Print[Unevaluated[AddToCompleteSquare[1]]," = ", AddToCompleteSquare[1]];

	   (* Complete the square in the denominator *)
	   Print["Denominator with completed square: "];
	   Denom = FeynDenom[1] /. {pr -> pr + AddToCompleteSquare[1]} // DiracAlgebra // Factor;
	   FeynDenom[1] = -( Denom  - Scal[pr, pr]) //Factor;
	   Print[Unevaluated[FeynDenom[1]],"   = ", FeynDenom[1]];
	   Return[FeynCoeff[0] (Scal[pr,pr] - FeynDenom[0])^(-FeynPow[0])];

	 ];






FeynParamAndMore2[{a___}, pr_, number_?NumberQ] :=     
    Module[{arg, PowSum, GammaProd, denom, fpPrefac, i, res, Denom},   

	  (*  a of  the  form  {den1, pow1, fp1}, {den2, pow2, fp2},  ... , {denN, powN, fpN}  *)
	  (*  Meaning: 1/den^pow1 1/den2^pow2*...*1/denN^powN, fpI are  the desired Feynman Variables *)
	  (*  PowSum :   Sum  of  the  powers; need  that  for  Gamma[ pow1 +   ... ]  *)
	  (*  GammaProd :   the  denominator  of  the  prefactor  *)
	  (*  denom  collects  the  final  denominator  *)
	  (*  fpPrefac  collects  the  prefactors  involving  Feynamn Parameters  *)

	  arg   =   {a};  
	  For[PowSum = 0; GammaProd = 1; fpPrefac = 1; denom = 0;    
	      i = 1, i <= Length[arg],  i++,   
	      PowSum    += arg[[i, 2]];   
	      GammaProd *= Gamma[arg[[i, 2]]];   
	      fpPrefac  *= arg[[i, 3]]^(arg[[i, 2]] - 1); 
	      denom     += arg[[i, 1]] * arg[[i, 3]];
	     ]; 

	  (*   res = Gamma[PowSum]/GammaProd fpPrefac denom^(-PowSum); *)   
	      
	  FeynCoeff[number, 1] = Gamma[PowSum]/GammaProd * fpPrefac;   
	  FeynDenom[number, 1] = denom;
	  FeynPow[number, 1]   = PowSum;  
       	  Print[Unevaluated[FeynCoeff[number, 1]]," = ", FeynCoeff[number, 1]];
	  Print[Unevaluated[FeynDenom[number, 1]]," = ", FeynDenom[number, 1]];
	  Print[Unevaluated[FeynPow[number, 1]],"   = ", FeynPow[number, 1]];

	   (* And Complete the Square *)
	   AddToCompleteSquare[number, 1] = CompleteTheSquare[FeynDenom[number, 1], pr];
	   Print[Unevaluated[AddToCompleteSquare[number, 1]]," = ", AddToCompleteSquare[number, 1]];

	   (* Complete the square in the denominator *)
	   Print["Denominator with completed square: "];
	   Denom = FeynDenom[number, 1] /. {pr -> pr + AddToCompleteSquare[number, 1]} // DiracAlgebra // Factor;
	   FeynDenom[number, 1] = -( Denom  - Scal[pr, pr]) //Factor;
	   Print[Unevaluated[FeynDenom[number, 1]],"   = ", FeynDenom[number, 1]];
	   Return[FeynCoeff[number,0] (Scal[pr,pr] - FeynDenom[number, 0])^(-FeynPow[number, 0])];

	 ];


(************************************************End of Tensor Reduction***********************************)

(***************************************** Gamma- and Betafunctions *****************************************)

(* Simplify expressions with Gamma functions *)
GammaExpand={ Gamma[a__] :> Gamma[Expand[a]]};
GammaSimplifyList = {	
    Gamma[n_?NumberQ + a_. ] :> (n-1+a) Gamma[n-1+a] /; n>=1,
    Gamma[n_?NumberQ - a_. ] :> (n-1-a) Gamma[n-1-a] /; n>=1,
    Gamma[n_?NumberQ + a_. ] :> 1/(n+a) Gamma[n+1+a] /; n<0,
    Gamma[n_?NumberQ - a_. ] :> 1/(n-a) Gamma[n+1-a] /; n<0};

GammaSimplify[expr_] := Module[{res},
			  res = expr/.GammaExpand;
			  res = res //.GammaSimplifyList;
			  Return[res];
			  ];


BetaToGamma = { Beta[a_, b_] :> (Gamma[a] Gamma[b])/Gamma[a + b]  };


(* Replace -1+x and 1-x by +/- OneMinus[x] may use:  epxr/.ReplOneMinus[x,y,u,v]         *)
ReplOneMinus[u__, v_]:= Union[ReplOneMinus[u], ReplOneMinus[v]];
ReplOneMinus[u_]:= {1-u :> OneMinus[u], -1 + u :> - OneMinus[u], (OneMinus[u] u)^p_. :> OneMinus[u]^p u^p};
ReplOneMinus[{u_,v__}] := Union[ReplOneMinus[{u}],ReplOneMinus[{v}]];
ReplOneMinus[{u_}]:= {1-u :> OneMinus[u], -1 + u :> - OneMinus[u], (OneMinus[u] u)^p_. :>OneMinus[u]^p u^p};



(* Beta integrals and remains *)
(* !!!! Make sure to always use expr//.BetaIntegrals (ReplaceRepeated!) !!!  *)
BetaIntegrals[u__, v_] := Union[BetaIntegrals[u], BetaIntegrals[v]];
BetaIntegrals[{u__, v_}] := Union[BetaIntegrals[{u}], BetaIntegrals[{v}]];
BetaIntegrals[{u_}] := {u^p_.  OneMinus[u]^q_.  :>  Beta[p + 1, q + 1]};
BetaIntegrals[u_] := {u^p_.  OneMinus[u]^q_.  :>  Beta[p + 1, q + 1]};

BetaRemains[u__, v_] := Union[BetaRemains[u], BetaRemains[v]];
BetaRemains[u_] := {u^p_.  :> Beta[p+1,1], OneMinus[u]^q_.  :> Beta[1, q + 1] };
BetaRemains[{u__, v_}] := Union[BetaRemains[{u}], BetaRemains[{v}]];
BetaRemains[{u_}] := {u^p_.  :> Beta[p+1,1], OneMinus[u]^q_.  :> Beta[1, q + 1] };


                
 (********************************************************
        
        Partial Fraction
 ****** *******************************************************)   
                
(*************************************************************

                     Deleting of Propagators

**************************************************************)

    deletprop:=
  { 


    AD[den[a___], m___, den[a___], den[b_, 0], m2___, den[b_, 0]] :> 0, 
    AD[den[a_, 0], m___, den[a_, 0], den[b___], m2___, den[b___]] :> 0,
    
    AD[den[a_, 0], m___, den[a_, 0], den[b___]] :> 0,
    AD[den[b___], m2___, den[b___], den[a_, 0]] :> 0,
    
    AD[den[a_, 0], den[b___], m2___, den[b___]] :> 0,
    AD[den[a___], den[b_, 0], m___, den[b_, 0]] :> 0,
    
    AD[den[a_, 0], den[b___]] :> 0,
    AD[den[b___], den[a_, 0]] :> 0,
    
    AD[den[a_, 0], m___, den[a_, 0]] :> 0,
    AD[den[a_, 0]] :> 0
    

}


delet = Dispatch[deletprop];
	    Deletprop=Dispatch[deletprop];

partq1 := 
    AD[l___, den[q1_, m1_], m___, den[q1_, m2_], r___] :> 
      AD[l, den[q1, m1],m, r]/(m1^2 - m2^2) - 
          AD[l,m, den[q1, m2],r]/(m1^2 - m2^2) /; FreeQ[m1, m2];

partq2:=Scal[q1_,q1_]*AD[l___,den[q1_,m1_],m___,den[q1_,m2_],r___]:>
1/(m1^2-m2^2)*m1^2*
AD[l,den[q1,m1],m,r]-
1/(m1^2-m2^2)*m2^2*
AD[l,m,den[q1,m2],r]/;FreeQ[m1,m2];

partq40:= 
AD[l___,den[q1_, m2_],m___,den[q2_, m3_],n___,den[q1_ + q2_, m1_], r___]* Scal[q1_, q2_]^k_. :> 
1/2*Scal[q1,q2]^(k-1)*AD[l, den[q1, m2],m, den[q2, m3],n, r] - 
1/2*Scal[q1,q2]^(k-1) AD[l, m, den[q2, m3],n, den[q1 + q2, m1],r] - 
1/2*Scal[q1,q2]^(k-1) AD[l, den[q1, m2],m,n,den[q1 + q2, m1],r] + 
1/2*Scal[q1,q2]^(k-1) AD[l, den[q1, m2],m,den[q2, m3],n,den[q1 + q2, m1], r](m1^2 - m2^2 - m3^2);


partq5:=Scal[q1_,q1_]^n_. AD[l___,den[q1_,m1_],r___]:>
Scal[q1,q1]^(n-1) AD[l,r]+Scal[q1,q1]^(n-1)*m1^2*AD[l,den[q1,m1],r];



partq6:=Scal[q1_,q2_] AD[l___,den[q1_+q2_,m1_],r___]:>
  1/2(Scal[q1+q2,q1+q2]-Scal[q1,q1]-Scal[q2,q2])*AD[l,den[q1+q2,m1],r];


PartialFractionOne[expr__] := 
  Module[{part, part1, part2, part3,xx,pippo},

	Clear[pippo,pippo2];
 part = expr; 

While
[
pippo=!=part, 
	pippo=part;
part = part //.{partq40,partq5,partq2,partq1,partq6};

part=Expand[part];
];
part=part//.AD[]:>0;

part=part//.AD[a___]:>Sort[AD[a]];

part=part//.AD[den[q1_, 0], m___, den[q2_, 0], m2___, den[q1_ + q2_, 0]] :> 0;
 Return[part];]




PartialFractionTwo[expr__] := 
  Module[{part, part1, part2, part3,xx,pippo},

	Clear[pippo,pippo2];
 part = expr; 

While
[
pippo=!=part, 
	pippo=part;
part = part //.{partq40,partq5,partq2,partq1,partq6};

part=Expand[part];
];
part=part//.AD[]:>0;
part=part//.AD[a_]:>0;
part=part//.AD[a___]:>Sort[AD[a]];
part=part//.delet;
part=part//.AD[den[q1, 0], m___, den[q2, 0], m2___, den[q1 + q2, 0]] :> 0;
 Return[part];]


   (**************************************End of partial fraction*******************************)

 
 
 
 (***************************Substitutions**********************************)

(********Substitutions for three masses, two different******************************************)




subax :={ (c_.)*AD[l___, den[(q1_)?LoopMomentumQ - 
       (q2_)?LoopMomentumQ, M_], r___] :> 
       ( c*AD[l, den[q1 + q2, M], r]/;FreeQ[c,q2])
//DiracAlgebra//Expand}

   


     subaxino:= G[i[m2_  ,n2_], i[m1_, n1_], i[m1_, n3_]] :> 
      G[i[m1, n1], i[m2, n2], i[m1,n3]]/;FreeQ[m1,m2]

   

(***************************Substitutions Rest**********************************)

Substitution[expr__] := 
  Module[{part},
part = expr; 
 part=part/.c_.*AD[den[(q1_)?LoopMomentumQ, m1_], m___, den[(q1_)?LoopMomentumQ, m1_], 
den[(q1_)?LoopMomentumQ + (q2_)?LoopMomentumQ, m2_], r___]
 :>
(c*AD[den[q1, m1], m, den[q1, m1], den[q1 + q2, m2], r]
/.q2 -> (q2 - q1) // DiracAlgebra // Expand);

 part=part/.c_.*AD[den[(q1_)?LoopMomentumQ, m1_], den[(q1_)?LoopMomentumQ + (q2_)?LoopMomentumQ, m2_], r___]
:>
(c*AD[den[q1, m1], den[q1 + q2, m2], r]
/.q2 -> (q2 - q1) // DiracAlgebra // Expand);

     part=part/.subax;



 part=part//.AD[a___]:>Sort[AD[a]];

 Return[part];]



  
(***************************End of Substitutions **********************************)

  


(**************************************************************************
              Simplifaction of Propagator structure
**************************************************************************)
     
forwardpower := {AD[a__, b___] :> AD[a*b], AD[a__, b___] -> AD[a*b]}

middlepower := AD[a___*b___] :> AD[a, b];

backpower := (AD[c___, Power[den[a__], n_], d___] :> 
AD[c, Power[den[a], n - 1], den[a], d] /; n > 1)
			     (*hier muss momentumq masse und impuls abgefragt werden*)
     power := {AD[l___, den[q1_, m1_]^(n_.), r___] :> AD[l,i[m1, n],r]}


sub1 := G[i[0, a1_], i[m2_, a2_], i[m3_, a3_]] :> 
G[i[m2, a2], i[m3, a3], i[0, a1]];
sub2 := G[i[m2_, a2_], i[0, a1_], i[m3_, a3_]] :> 
G[i[m2, a2], i[m3, a3], i[0, a1]];
sub3 := 
G[i[Global`MT, a2_], i[MW_, a3_], i[0, a1_]] :> G[i[MW, a3], i[Global`MT, a2], i[0, a1]];
sub5 := G[i[Sqrt[x1_] M1_, a2_], i[M1_, a3_], i[0, a1_]] :> 
G[i[M1, a3], i[Sqrt[x1]M1, a2], i[0, a1]];
sub6 := G[i[M2_, a2_], i[Global`M1, a3_], i[0, a1_]] :> 
G[i[Global`M1, a3], i[M2, a2], i[0, a1]];
charmrule := 
G[i[0, a1_], i[0, a2_], i[m3_, a3_]] :> G[i[m3, a3], i[0, a2], i[0, a1]]
subaxino:= G[i[m2_  ,n2_], i[m1_, n1_], i[m1_, n3_]] :> 
      G[i[m1, n1], i[m2, n2], i[m1,n3]]/;FreeQ[m1,m2]
      
subrule={charmrule,sub1,sub2,sub3,sub5,sub6,subaxino};

SimplifyPropagator[expr__] := 
  Module[{part, part1, part2, part3,xx,pippo},
	part = expr; 
	part=part //. AD[a___] :> Sort[AD[a]]; 
	part=part//.delet;
	part = part //. forwardpower; 
	part = part //. middlepower; 
	part = part //. power;
	part=part/. AD[a_, b_, c_] :> G[a, b, c];
	part=part//.subrule;
	Return[part];]
			      

  
  
  (**************************************************************************
              End of Simplifaction of Propagator structure
**************************************************************************)
  






  

(*****************************************************************************

               Calculation of Two_Loop-Diagrams

 *****************************************************************************)







(*****************************Recurrence relations***************************)





				 

recu1 := C[i[m1_, a1_], i[m2_, a2_], i[0, a3_]] :> 
    1/((-1 + a1)  (1 - Global`x1))((1 - a1 - a2 - a3 + 4 - 
                  2 *eps + (-1 + a1 - a3) *Global`x1) C[i[m1, -1 + a1], i[m2, a2], 
                i[0, a3]] + 
            Global`x1*a2  (C[i[m1, -2 + a1], i[m2, 1 + a2], i[0, a3]] - 
                  C[i[m1, -1 + a1], i[m2, 1 + a2], i[0, -1 + a3]])) /; (a1 > 
            1 && a1 > a2 && a1 > a3 && a2 > 0)




recu2 := C[i[m1_, a1_], i[m2_, a2_], i[0, a3_]] :>
 -1/Global`x1*((-1 + a2) (1 - Global`x1))*( ((-1 + a2 - a3 + (1 - a1 - a2 - a3 + 4 - 2 eps) Global`x1) 
C[i[m1, a1], i[m2, -1 + a2], i[0, a3]] + 
              a1 (C[i[m1, 1 + a1], i[m2, -2 + a2], i[0, a3]] - 
                    C[i[m1, 1 + a1], i[m2, -1 + a2], i[0, -1 + a3]]))) /; 
      a2 > 1 && a2 >= a1 && a2 >= a3 && a1 > 0



recu3 := C[i[m1_, a1_], i[m2_, a2_], i[0, a3_]] :> 
    1/((a3 - 1)*(1 - Global`x1)^2)*
(((Global`x1 + 1)(-4 + 2 eps) + 2*a2 + (1+3 Global`x1)*(a3 - 1)) C[i[m1, a1], i[m2, a2], i[0, a3 - 1]] + 2 *Global`x1* a2 (C[i[m1, a1], i[m2, a2 + 1], i[0, a3 - 2]] - 
            C[i[m1, a1 - 1], i[m2, a2 + 1], i[0, a3 - 1]]) +
(a3 - 1)(1 - Global`x1) (C[i[m1, a1], i[m2, a2 - 1], i[0, a3]] - 
                  C[i[m1, a1 - 1], i[m2, a2], i[0, a3]])) /; (a3 > 1 && 
          a3 >= a1 && a3 >= a2 && a1 > 0 && a2 > 0)






recrulealt:=
G[i[m1_, 1], i[m2_, 1], i[0, 1]] :> 
  Normal[Series[(Pi^4*m1^2*N2[m1]*((-1 - Global`x1)/eps^2 + 
       (2*Global`x1*log[Global`x1])/eps + (1 - 2*Global`x1)*log[Global`x1]^2 + 
       2*(1 - Global`x1)*PoLi2[1 - 1/Global`x1]))/(2*(1 - 2*eps)*
      (1 - eps)), {eps, 0, 0}]]

(*first equation of (59) Misiak*)

rec11 := G[i[m1_, a1_], i[m2_, a2_], i[0, a3_]] :> 
    1/(m1^2*(-1 + a1) (1 - Global`x1))((1 - a1 - a2 - a3 + 4 - 
                  2 eps + (-1 + a1 - a3) Global`x1) G[i[m1, -1 + a1], i[m2, a2], 
                i[0, a3]] + 
            Global`x1 a2 (G[i[m1, -2 + a1], i[m2, 1 + a2], i[0, a3]] - 
                  G[i[m1, -1 + a1], i[m2, 1 + a2], i[0, -1 + a3]])) /; 
      a1 > 1 && a1 > a2 && a1 > a3 && a2 > 0



rec22 := G[i[m1_, a1_], i[m2_, a2_], 
      i[0, a3_]] :> -1/((m1^2*
                Global`x1*(-1 + a2) (1 - Global`x1)))(((-1 + a2 - 
                    a3 + (1 - a1 - a2 - a3 + 4 - 2 eps) Global`x1) G[i[m1, a1], 
                  i[m2, -1 + a2], i[0, a3]] + 
              a1 (G[i[m1, 1 + a1], i[m2, -2 + a2], i[0, a3]] - 
                    G[i[m1, 1 + a1], i[m2, -1 + a2], i[0, -1 + a3]]))) /; 
      a2 > 1 && a2 >= a1 && a2 >= a3 && a1 > 0

      rec33 := G[i[m1_, a1_], i[m2_, a2_], i[0, a3_]]:> 1/(m1^2(a3-1)*(1-Global`x1)^2)*(((Global`x1 + 1) (-4 + 2 eps) + 2 a2 + (1 + 3 Global`x1) (a3 - 1)) G[i[m1, a1], i[m2, a2], 
          i[0, a3 - 1]] + 
      2 Global`x1 a2 (G[i[m1, a1], i[m2, a2 + 1], i[0, a3 - 2]] - 
            G[i[m1, a1 - 1], i[m2, a2 + 1], i[0, a3 - 1]]) + (a3 - 1) (1 - 
            Global`x1) (G[i[m1, a1], i[m2, a2 - 1], i[0, a3]] - 
            G[i[m1, a1 - 1], i[m2, a2], i[0, a3]])) /; 
  a3 > 1 && a3 >= a1 && a3 >= a2 && a1 > 0 && a2 > 0

(*to change the notation
m1^2-m2^2= (1-Global`x1)m1^2
m1^2+m2^2=()m1^2
 *)
delta[M1_,M2_,M3_]=2(M1^2 M2^2 +M1^2 M3^2 +M2^2 M3^2)-(M1^4 +M2^4+M3^4)

recu11:=G[i[m1_, a1_], i[m2_, a2_], i[m3_, a3_]] :>
1/((a1-1) m1^2*delta[m1,m2,m3])((a2 (m1^2 -m3^2)(m1^2-m2^2+m3^2)
				 +a3(m1^2 -m2^2)(m1^2+m2^2-m3^2)+d m1^2 (-m1^2+m2^2+m3^2)-(a1-1) delta[m1,m2,m3]) G[i[m1,a1-1],i[m2,a2],i[m3,a3]]
				+a2 m2^2 (m1^2-m2^2+m3^2)
				(G[i[m1,a1-1],i[m2,a2+1],i[m3,a3-1]]-G[i[m1,a1-2],i[m2,a2+1],i[m3,a3]])
+a3 m3^2 (m1^2+m2^2-m3^2)
(G[i[m1,a1-1],i[m2,a2-1],i[m3,a3+1]]-G[i[m1,a1-2],i[m2,a2],i[m3,a3+1]]))/; a1 > 1 && a1 > a2 && a1 > a3 && a2 > 0

(*recu11 corresponds to rec11 in the limit m3->0*)
(*recu22 corresponds to rec22 in the limit m3->0*)
(*recu33 corresponds to rec22 in the limit m3->0*)
recu22:=G[i[m1_, a1_], i[m2_, a2_], i[m3_, a3_]] :>
1/((a2-1) m2^2*delta[m1,m2,m3])((a1 (m2^2 -m3^2)(m2^2-m1^2+m3^2)
				 +a3(m2^2 -m1^2)(m1^2+m2^2-m3^2)+d m2^2 (-m2^2+m1^2+m3^2)-(a2-1) delta[m1,m2,m3]) G[i[m1,a1],i[m2,a2-1],i[m3,a3]]
				+a1 m1^2 (m2^2-m1^2+m3^2)
				(G[i[m1,a1+1],i[m2,a2-1],i[m3,a3-1]]-G[i[m1,a1+1],i[m2,a2-2],i[m3,a3]])
+a3 m3^2 (m1^2+m2^2-m3^2)
(G[i[m1,a1-1],i[m2,a2-1],i[m3,a3+1]]-G[i[m1,a1],i[m2,a2-2],i[m3,a3+1]]))/; 
      a2 > 1 && a2 >= a1 && a2 >= a3 && a1 > 0

recu33:=G[i[m1_, a1_], i[m2_, a2_], i[m3_, a3_]] :>
1/((a3-1) m3^2*delta[m1,m2,m3])((a2 (m3^2 -m1^2)(m1^2-m2^2+m3^2)
				 +a1(m3^2 -m2^2)(m3^2+m2^2-m1^2)+d m3^2 (-m3^2+m2^2+m1^2)-(a3-1) delta[m1,m2,m3]) G[i[m1,a1],i[m2,a2],i[m3,a3-1]]
				+a2 m2^2 (m1^2-m2^2+m3^2)
				(G[i[m1,a1-1],i[m2,a2+1],i[m3,a3-1]]-G[i[m1,a1],i[m2,a2+1],i[m3,a3-2]])
+a1 m1^2 (m3^2+m2^2-m1^2)
(G[i[m1,a1+1],i[m2,a2-1],i[m3,a3-1]]-G[i[m1,a1+1],i[m2,a2],i[m3,a3-2]]))/; 
  a3 > 1 && a3 >= a1 && a3 >= a2 && a1 > 0 && a2 > 0


recu333:=G[i[m1_, a1_], i[m2_, a2_], i[m3_, a3_]] :>
1/((a3-1) m3^2*delta[m1,m2,m3])((a1 (m3^2 -m2^2)(m2^2-m1^2+m3^2)
				 +a2(m3^2 -m1^2)(m1^2+m3^2-m2^2)+d m3^2 (-m3^2+m1^2+m2^2)-(a3-1) delta[m1,m2,m3]) G[i[m1,a1],i[m2,a2],i[m3,a3-1]]
				+a1 m1^2 (m2^2-m1^2+m3^2)
				(G[i[m1,a1+1],i[m2,a2-1],i[m3,a3-1]]-G[i[m1,a1+1],i[m2,a2],i[m3,a3-2]])
+a2 m2^2 (m1^2+m3^2-m2^2)
(G[i[m1,a1-1],i[m2,a2+1],i[m3,a3-1]]-G[i[m1,a1],i[m2,a2+1],i[m3,a3-2]]))/; a3 > 1 && a3 >= a1 && a3 >= a2 && a1 > 0 && a2 > 0

reconaxino := {G[i[m1_, a1_], i[m2_, a2_], i[0, a3_]] :> 
      0 /; a1 <= 0 || a2 <= 0, 
    G[i[m1_, a1_], i[m2_, a2_], i[0, 0]] :> AD[i[m1, a1], i[m2, a2]],
 G[l___, i[m_, 0], r___] :> AD[l, r]}




(*********************************End of recurrence relations*************************************************************************)





(**************** Two Loop Integration****************************)









(*Two Loop with Two Masses, no propagator without mass*)

(*all formula: Davydychev And Tausk, 1992*)
    
       (*Davydychev And Tausk: eq. 3.14*)
       (*z=Global`x1/4*)

       (* if Global`x1<<1 we have*)
         phix1=
         4 Sqrt[Global`x1/(4-Global`x1)] *ClausenCl2[2 ArcSin[Sqrt[Global`x1]/2]]
        
         lambda[x_]:=Sqrt[1-4/x]
         (* if Global`x1>1 we have*)
           
         phix2=1/lambda[Global`x1]        * (
         -4 * PoLi2[(1-lambda[Global`x1])/2]
         + 2 *(log[(1-lambda[Global`x1])/2])^2- (log[Global`x1])^2
         +Pi^2/3
         )
                        
                        
ca1 = -((1 + Global`x1/2)/eps^2) + (Global`x1*log[Global`x1])/eps - 
(1/2)*Global`x1*log[Global`x1]^2 + (2 - Global`x1/2) phix1;


ca2 = -((1 + Global`x1/2)/eps^2) + (Global`x1*log[Global`x1])/eps - 
(1/2)*Global`x1*log[Global`x1]^2 + (2 - Global`x1/2) phix2;


caa1=ca1/((1-eps)(1-2eps));

caa2=ca2/((1-eps)(1-2eps));



(*Davydev und Tausk (3.15*)
scaa=(-(3/(2*eps^2)) + 2*Sqrt[3]*ClausenCl2[Pi/3])/((1 - 2*eps)*(1 - eps))


ClausenCl[n_Integer?OddQ,x_]:=(PolyLog[n, E^(-I x)] + PolyLog[n, E^(I x)])/2
ClausenCl[n_Integer?EvenQ,x_] :=(PolyLog[n, E^(-I x)] - PolyLog[n, E^(I x)])I/2


clausen:=G[i[m1_, 1], i[m2_, 1], i[m1_, 1]] :> 
If[
Global`x<1,
Normal[Series[Pi^4*m1^2*N2[m1]*caa1,{eps,0,0}]],
 Normal[Series[Pi^4*m1^2*N2[m1]*caa2,{eps,0,0}]],
 Print["Declare the mass hierarchy of ", m1, " and ", m2, " by tying in an explicit value for x (=",m2^2/m1^2,")." ] 
    Abort[];
]





simpleclausen:=G[i[m1_, 1], i[m1_, 1], i[m1_, 1]] :> 
Normal[Series[Pi^4*m1^2*N2[m1]*scaa,{eps,0,0}]]


ScalIntTwoThreeMasses[expr__] := 
    Module[{part, part1, part2, part3, xx, pippo}, Clear[a,b,pippo, pippo2];
      part = expr;
           part=part/.Global`MT->Global`M2;
   While[pippo =!= part, 
        pippo = part;
        part = part /.recurrence;
        part=part/.dToEps;
        part = part /. reconaxino;];
part=part//.simpleclausen;   
part=part//.subaxino; 
part=part//.clausen;
part=part//.facrule; 
  part=part/.nerules;  
        part = Normal[Series[part, {eps, 0, 0}]];
      
      Return[part];]
      
ScalIntTwo[expr__] := 
Module[{part, part1, part2, part3, xx, pippo}, Clear[pippo, pippo2];
part = expr;
part=part/.AD[]->0;
part = part /.subrule;
part=part/.twoequalmasses;
part = part/.twozeromasses;
While[pippo =!= part, 
pippo = part;
part = part /. recurrenceb;
part = part /. reconalt;];
part = part /. recrulealt;
part = part /. facrule;
 part = part // Expand;
part = Normal[Series[part, {eps, 0, 0}]];
part=part/.nerules;
Return[part];]

ScalIntTwoFull[expr__] := 
Module[{part, part1, part2, part3, xx, pippo}, Clear[pippo, pippo2];
part = expr;
part=part//.AD[]->0;
part=part//.Global`constantsM1M2;
part=part//.Global`MT->Global`M2;
part = part //. subrule;
While[pippo =!= part, pippo = part;
     part = part /. recurrenceb;
part = part /. reconalt;];
part = part /. recrulealt;
part = part /. facrule;
 part = part // Expand;
part = Normal[Series[part, {eps, 0, 0}]];
part=part/.nerules;
Return[part];]




    
(*change all factorising two-loop-integrals to this form*)
     
recurrence:=Dispatch[{recu11, recu22, recu33}]
recurrenceb:=Dispatch[{rec11, rec22, rec33}]



Recurrence[expr__] := 
Module[{part, part1, part2, part3, xx, pippo}, Clear[pippo, pippo2];
part = expr;
part=part//.AD[]->0;
part=part//.Global`constantsM1M2;
part=part//.Global`MT->Global`M2;
part = part //. subrule;
While[pippo =!= part, pippo = part;
part=part/.reconaxino;
part=part/.recurrence;
part = part /. recurrenceb;
part = part /. reconalt;];
Return[part];]


facruleneu := AD[i[M1_, a_], i[M2_, b_]] :> 
-(Pi^4*(-1)^a*(-1)^b*Pochhammer[1 + eps, a - 3]*Pochhammer[1 + eps, b - 3]*
(1 - log[Global`x1]*eps + (1/2)*log[Global`x1]^2*eps^2)*N2[M1]*(M1^2)^(2 - a)*
(Global`x1*M1^2)^(2 - b))/((a - 1)!*(b - 1)!)
     
     facrule := 
   AD[i[m1_, a_], i[m2_, b_]] :> -(Pi^4*(-1)^a*(-1)^b*Pochhammer[1 + eps, a - 3]*
Pochhammer[1 + eps, b - 3]*Ne[m1]*Ne[m2]*(m1^2)^(2 - a)*(m2^2)^(2 - b))/
((a - 1)!*(b - 1)!)
     
     
     reconalt := {G[i[m1_, a1_], i[m2_, a2_], i[0, a3_]] :> 
0 /; a1 <= 0 || a2 <= 0, 
		    G[i[m1_, a1_], i[m2_, a2_], i[0, 0]] :> AD[i[m1, a1], i[m2, a2]]}

(*Ne[m1]=Pi^2*(1- eps kappa+O(eps)^2*)
     
     (*AD[i[m1,a]]integral d^D q 1/((q^2-m1^2)^a), no renormalization factor 1/(2 Pi^4 is included*)
 
     
     
    
     facruleone:= AD[i[m1_, a_]] :> Pi^2*I*(((-1)^a*Pochhammer[1 + eps, a - 3]*Ne[m1]*(m1^2)^(2 - a))/(a - 1)!)
     
     
     
 (**************Rules for simplifying prefactor structure*******************************)  
     
(*General Rules*)
nerule := Ne[M1_]*Ne[Sqrt[Global`x1]*M1_] :> (1 - log[Global`x1]*eps + (1/2)*log[Global`x1]^2*eps^2)*N2[M1]
nerule2 := Ne[Sqrt[Global`x1]*M1_]*Ne[Sqrt[Global`x1]*M1_] :> (1 - 2*eps*log[Global`x1] + 2*eps^2*log[Global`x1]^2)*     N2[M1]
nerule22 := Ne[Sqrt[Global`x1]*M1_]^2 :> (1 - 2*eps*log[Global`x1] + 2*eps^2*log[Global`x1]^2)*N2[M1]
nerule3 := N2[Sqrt[Global`x1]*M1_] :> (1 - 2*eps*log[Global`x1] + 2*eps^2*log[Global`x1]^2)*N2[M1]
nerule4 := Ne[M1_]^2 :> N2[M1]  


(*Rules for Global Masses M1 and M2*)

nerule5:=  Ne[Global`M1]* Ne[Global`M2] :> (1 -log[Global`x1]*eps +  (1/2)*log[Global`x1]^2*eps^2)*N2[Global`M1]
nerule6:=N2[Global`M2] :> (1 - 2*eps*log[Global`x1] + 2*eps^2*log[Global`x1]^2)*N2[Global`M1]

nerules=Dispatch[{
        Ne[M1_]*Ne[Sqrt[Global`x1]*M1_] :> (1 - log[Global`x1]*eps + (1/2)*log[Global`x1]^2*eps^2)*N2[M1],  
        Ne[Sqrt[Global`x1]*M1_]*Ne[Sqrt[Global`x1]*M1_] :> (1 - 2*eps*log[Global`x1] + 2*eps^2*log[Global`x1]^2)*N2[M1],
        N2[Sqrt[Global`x1]*M1_] :> (1 - 2*eps*log[Global`x1] + 2*eps^2*log[Global`x1]^2)*N2[M1],
        Ne[M1_]^2 :> N2[M1],
        Ne[Global`M1]* Ne[Global`M2] :> (1 -log[Global`x1]*eps +  (1/2)*log[Global`x1]^2*eps^2)*N2[Global`M1],
        N2[Global`M2] :> (1 - 2*eps*log[Global`x1] + 2*eps^2*log[Global`x1]^2)*N2[Global`M1]
        }]





fac=Dispatch{facrule,nerule,nerule2,nerule22,nerule4}






 (**************End of Rules for simplifying prefactor structure*******************************)  
     
     
(*******************End of Two Loop integration**************************************************)     

     
     
     
     
     
  (**************Scalar One LoopIntegral*******************************)  
         
     
     
ScalIntOne[expr__] := 
Module[{part, part1, part2, part3, xx, pippo}, Clear[a,pippo, pippo2];
part = expr;
part = part //. facruleone;
part = part // Expand;
part = part // Cancel;       
part = Normal[Series[part, {eps, 0, 1}]];

       Return[part];]




(*Two Loop Integral with two zero masses*)

twozeromasses:=G[i[m1_, a1_],i[0, a2_], i[0, a3_]] :> 
Normal[Series[(Pi^4*N2[m1]*(m1^2)^(4 - a1 - a2 - a3)*(-1)^(1 + a1 + a2 + a3)*
Pochhammer[1 + 2*eps, a1 + a2 + a3 - 5]*Pochhammer[1 + eps, a2 + a3 - 3]*
Pochhammer[1 - eps, 1 - a2]*Pochhammer[1 - eps, 1 - a3])/
((-1 + a1)!*(-1 + a2)!*(-1 + a3)!*(1 - eps)*(1 - (1/3)*Pi^2*eps^2)), {eps, 0, 0}]];


TwoZero[expr__] := 
      Module[{part}, Clear[a,pippo, pippo2];
part = expr;
part = part /.twozeromasses;
part=part/.nerule3;
part = part // Expand;
part = part // Cancel;
(* part = Normal[Series[part, {eps, 0, 0}]];*)

     Return[part];]

(*Two Loop Integral with two equal masses*)
      
      
      
twoequalmasses:=G[i[m1_,a1_],i[m1_,a2_],i[0,a3_]]:>Normal[Series[Pi^4*(N2[m1] (m1^2)^(4-a1-a2-a3) (-1)^(1+a1+a2+a3) Pochhammer[2-eps,-a3] Pochhammer[1+eps,-3+a1+a3] Pochhammer[1+eps,-3+a2+a3])/((-1+a1)! (-1+a2)! Pochhammer[a1+a2+a3-4+2 eps,a3]),{eps,0,0}]]



TwoEqual[expr__] := 
Module[{part, part1, part2, part3, xx, pippo}, Clear[a,pippo, pippo2];
part = expr;
part = part /.twoequalmasses;
part = part //.nerule3;
      part = part // Expand;
part = part // Cancel;
(* part = Normal[Series[part, {eps, 0, 0}]];*)
     
     
     
     Return[part];]

















man1 := Dirac[l___, q1_, q1_, r___] :> q1 \:05e0q1*Dirac[l, r] /; MomentumQ[q1]

man2 := Dirac[ll___, m1_, m2_, m1_, rr___] :> 
(-(m1 \:05e0m1))*Dirac[ll, m2, rr] + 2*m1 \:05e0m2*Dirac[ll, m1, rr] /; 
MomentumQ[m1] && MomentumQ[m2]

man22 := Dirac[l___, q1_, m___, q2_, q1_, r___] :> 
2*q1 \:05e0q2*Dirac[l, q1, m, r] - Dirac[l, q1, m, q1, q2, r] /; 
FreeQ[(FreeQ[{m}, #1, Infinity] & ) /@ {q2, Z2, Z3}, False] && MomentumQ[q1] && 
MomentumQ[q2]

man3 := Dirac[l___, q2_, q1_, lor1_, q2_, q1_, r___] :> 
2*Dirac[l, q2, r]*lor1 \:05e0q2*q1 \:05e0q1 - 2*Dirac[l, lor1, q1, q2, r]*q1 \:05e0q2 - 
     4*Dirac[l, q1, r]*lor1 \:05e0q2*q1 \:05e0q2 + 4*Dirac[l, lor1, r]*(q1 \:05e0q2)^2 + 
2*Dirac[l, q1, r]*lor1 \:05e0q1*q2 \:05e0q2 - Dirac[l, lor1, r]*q1 \:05e0q1*q2 \:05e0q2 /; 
IndexQ[lor1] && MomentumQ[q1] && MomentumQ[q2]

man := {man3, man2, man1}; 
manneu := {man3, man2, man1, man22}; 


PartialFractionTwo[expr__] := 
  Module[{part, part1, part2, part3,xx,pippo},

	Clear[pippo,pippo2];
 part = expr; 

While
[
pippo=!=part, 
	pippo=part;
part = part //.{partq40,partq5,partq2,partq1};

part=Expand[part];
];
part=part//.AD[]:>0;
part=part//.AD[a_]:>0;
part=part//.AD[a___]:>Sort[AD[a]];
part=part//.delet;
part=part//.AD[den[q1, 0], m___, den[q2, 0], m2___, den[q1 + q2, 0]] :> 0;
 Return[part];]


finalpart:=
  a_.*Scal[q1_,q2_]^n_:>
    PartialFractionTwo[a*Scal[q1,q2]^n]/;LoopMomentumQ[q1]&&LoopMomentumQ[q2]


      FinalPart[expr_]:=Module[{part},part=expr/.finalpart;Return[part];]


denrule1:=den[-q1+k1_,m_]:>den[q1-k1,m]/;FreeQ[k1, q2]
denrule2:=den[-q2+k1_,m_]:>den[q2-k1,m]/;FreeQ[k1, q1]
denrule3:=den[-q1-q2+k1_,m_]:>den[q1+q2-k1,m]
denrule4:=den[-q1+q2+k1_,m_]:>den[q1-q2-k1,m]/;FreeQ[k1, q1]&&FreeQ[k1, q1]
denrule5:=den[-q1+q2,m_]:>den[q1-q2,m]



denrulesub:={
den[-q1___+k1_,m_]:>
den[q1-k1,m]/;FreeQ[k1, _?LoopMomentumQ] && LoopMomentumQ[q1],
den[-q1_+q2_+k_,m_]:>den[q1-q2-k,m]/;LoopMomentumQ[q1]&&LoopMomentumQ[q2],
den[-q1_+q2_,m_]:>den[q1-q2,m]/;LoopMomentumQ[q1]&&LoopMomentumQ[q2]
}



subden:={
den[-q1___+k1_,m_]:>
den[q1-k1,m]/;FreeQ[k1, _?LoopMomentumQ] && LoopMomentumQ[q1],
den[-q1_+q2_+k_,m_]:>den[q1-q2-k,m]/;LoopMomentumQ[q1]&&LoopMomentumQ[q2],
				       den[-q1_+q2_,m_]:>den[q1-q2,m]/;LoopMomentumQ[q1]&&LoopMomentumQ}



  
  
  
  
  
  
  
(******************************************** Misc. stuff ********************************************)

dToEps = {d -> 4 - 2 eps};

ScalToTimes = { Scal -> Times };

FourFeyn[list1_?ListQ, list2_?ListQ]:=
    Module[{ala=0,bla=0,cla=0},
      ala=Length[list1];
      bla=Length[list2];
      If[ala==4,
        If[bla==3,      
          FourFeynFactor[1]=6*list2[[2]]*list2[[3]]*list2[[3]];
          Print["FourFeynFactor[1] = ",FourFeynFactor[1]];
          FourDenom[1]=Power[list1[[1]]*list2[[1]]*list2[[2]]*list2[[3]]+list1[[2]]*(1-list2[[1]])*list2[[2]]*list2[[3]]+list1[[3]]*(1-list2[[2]])*list2[[3]]+list1[[4]]*(1-list2[[3]]),1];
          FourPower[1]=4;
          Print["FourDenom[1] = ", FourDenom[1]],
          Print["Not the right amount of Feynman parameters! Exiting..."]],
      Print["Not the right amount of denominators! Exiting..."];];
      ];


(* Convert g-strong to alpha-strong and vice versa *)

gsToAlphas = {gs^n_ :> 4 Pi \[Alpha]s gs^(n-2)  /;  n >=  2};

AlphasTogs = {\[Alpha]s :> gs^2/(4 Pi)};


   
     selfrule= Dirac[R, k1] -> (EL*(1 - (2*SW^2)/3)*Dirac[R, lor1])/(2*CW*SW)
     
     zlepton=+EL*CW^2/(2*SW*CW*MW^2)*((1/2-2*SW^2) *Dirac[lor1]-1/2 Dirac[lor1,gamma5]);
     
     

     
  
  
 Protect["Integrals`*"];


  End[];    (* "`Private`" *) 

EndPackage[];

(*************************************************************************)
     
