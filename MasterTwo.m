
Off[General::spell, General::spell1]
Off[Syntax::"stresc"]
Off[General::shdw]

env := Print["Date: ", Date[]\[LeftDoubleBracket]3\[RightDoubleBracket], ".",
       Date[]\[LeftDoubleBracket]2\[RightDoubleBracket], ".", 
      Date[]\[LeftDoubleBracket]1\[RightDoubleBracket], "     Time: ", 
      Date[]\[LeftDoubleBracket]4\[RightDoubleBracket], ":", 
    If[Date[][[5]]<10,StringJoin[ToString[0],ToString[Date[][[5]]]], Date[][[5]]]];

env

Print["Before usage you have to declare all occuring Lorentz indices, masses and momenta!"]
BeginPackage["MasterTwo`",{"Fermions`","Integrals`"}]
MasterTwoInfo[] := (
    Print["\n *** \t MasterTwo \t ***"];
    Print["MasterTwo is a package that can simplify Dirac expressions 
          in d dimensions with an anticommuting Gamma5. Furthermore 
          it can Taylor expand and tensor reduce the integrands as 
          well as integrate scalar two-loop integrals with up to two
          different mass scales."];
    Print["\n Available commands:"];
    Print[" (For any of these commands '?command' gives a short
             description of its functionaity.)\n"];
    Print[Names["Fermions`*"]];
    Print[Names["Integrals`*"]];  
		   );
EndPackage[];

