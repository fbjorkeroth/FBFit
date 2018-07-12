(* ::Package:: *)

BeginPackage["FBFit`CalculateParameters`",{"FBFit`MixingParameterTools`MPT3x3`"}];

(*FBCalculateParametersL::usage ="takes the neutrino mass matrix mnu and charged lepton Yukawa matrix Ye as input and calculates the mixing angles, phases, and masses.";
FBCalculateParametersQ::usage ="takes the quark Yukawa matrices Yu, Yd as input and calculates the mixing angles, phase, and masses.";*)

FBCalculateParameters::usage ="takes the neutrino mass matrix mnu and fermion Yukawa matrices Yu, Yd, and Ye as input and calculates the various mixing angles, phases, and masses.";
FBGetPulls::usage="takes as input the result of running FBCalculateParameters and generates a list containing the pulls (i.e. deviation from best fit).";
FBChiSq::usage="takes the sum of squares of pulls to give the final \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\) value.";


Begin["`Private`"];

FBCalculateParameters[Yu_?MPTNumericMatrixQ,Yd_?MPTNumericMatrixQ,mnu_?MPTNumericMatrixQ,Ye_?MPTNumericMatrixQ]:=Flatten@Join[
	calculateParametersQ[Yu,Yd],
	calculateParametersL[mnu,Ye]
];

FBCalculateParameters[mat1_?MPTNumericMatrixQ,mat2_?MPTNumericMatrixQ,type_?StringQ]:=Module[{f},
	f=Switch[type,
		"Q",calculateParametersQ,
		"L",calculateParametersL
	];
	f[mat1,mat2]
];

FBGetPulls[calcparam_?(VectorQ[#,NumericQ]&),bestfit_?(VectorQ[#,NumericQ]&),errors_?(VectorQ[#,NumericQ]&)]:=
	Table[(calcparam[[i]]-bestfit[[i]])/errors[[i]],{i,Length@calcparam}];

FBChiSq[pulls_?(VectorQ[#,NumericQ]&)]:=Sum[pulls[[i]]^2,{i,Length@pulls}];

calculateParametersL[mnu_?MPTNumericMatrixQ,Ye_?MPTNumericMatrixQ]:=Module[{paraml,output},
	paraml=MNSParameters[mnu,Ye];
	output=List[];
	AppendTo[output,paraml[[1,1;;4]]/Degree];
	AppendTo[output,paraml[[2,#]]^2-paraml[[2,1]]^2&/@{2,3}];
	AppendTo[output,paraml[[3]]];
	Flatten@output
];

calculateParametersQ[Yu_?MPTNumericMatrixQ,Yd_?MPTNumericMatrixQ]:=Module[{paramq,output},
	paramq=CKMParameters[Yu,Yd];
	output=List[];
	AppendTo[output,paramq[[1]]/Degree];
	AppendTo[output,paramq[[2]]];
	AppendTo[output,paramq[[3]]];
	Flatten@output
];

End[];
EndPackage[];
