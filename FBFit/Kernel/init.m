(* Wolfram Language Init File *)

Get[ "FBFit`MixingParameterTools`MPT3x3`"]
Get[ "FBFit`CalculateParameters`"]
Get[ "FBFit`BestFitsAndErrors`"]
Get[ "FBFit`SetOptions`"]
Get[ "FBFit`Analysis`"]
Get[ "FBFit`FBFit`"]

Print["Package FBFit` loaded.\nFor an overview of the package functions, see the ",
	Hyperlink["documentation","paclet:FBFit",BaseStyle->Bold],
	".\nSeveral variables must be defined manually, see the ",
	Hyperlink["tutorial","paclet:FBFit/tutorial/A simple MCMC fit",BaseStyle->Bold]," for details."
];
