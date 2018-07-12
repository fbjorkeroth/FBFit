(* ::Package:: *)

(* ::Title:: *)
(*Simple MCMC*)


Get["FBFit`"];

FBLoadModel["models/C1.m"];
nMCMC=10000;
FBSetOptions[
	"Model"->"SM",
	"SaveOutput"->False,
	"ThinningSaveFile"->1,
	"BurnIn"->0
];
FBLoadBestFitsAndErrors[];
\[Theta]0=FBSetSeed[1];
l=FBMonteCarlo[nMCMC,\[Theta]0]//AbsoluteTiming



