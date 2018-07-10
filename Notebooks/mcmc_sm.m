(* ::Title:: *)
(*Simple MCMC*)

AppendTo[$Path,ParentDirectory[NotebookDirectory[]]];
Get["FBFit`"];

FBLoadModel["models/model.m"];

nMCMC=100;

FBSetOptions[
	"Model"->"SM",
	"BurnIn"->0,
	"SaveOutput"->False,
	"ExcludeParameters"->Range[11,19],
	"ThinningSaveFile"->Ceiling[nMCMC/50000]
];

FBLoadBestFitsAndErrors[];

\[Theta]0=FBSetSeed[];

l=FBMonteCarlo[nMCMC,\[Theta]0];