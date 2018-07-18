(* ::Package:: *)

(* ::Title:: *)
(*Simple MCMC*)


Get["FBFit`"];

FBLoadModel["models/model.m"];
nMCMC=10000;
FBSetOptions[
	"Model"->"SM",
	"SaveOutput"->True,
	"ThinningSaveFile"->1,
	"BurnIn"->1000,
	"Sector"->"Q"
];
FBLoadBestFitsAndErrors[];
\[Theta]0=FBSetSeed[];
l=FBMonteCarlo[nMCMC,\[Theta]0];



Beep[]
