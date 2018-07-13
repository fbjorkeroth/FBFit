(* ::Package:: *)

(* ::Title:: *)
(*Simple MCMC*)


Get["FBFit`"];

FBLoadModel["models/C8.m"];
nMCMC=50000;
FBSetOptions[
	"Model"->"SM",
	"SaveOutput"->True,
	"ThinningSaveFile"->1,
	"BurnIn"->10000,
	"Sector"->"Q"
];
FBLoadBestFitsAndErrors[];
\[Theta]0=FBSetSeed[];
l=FBMonteCarlo[nMCMC,\[Theta]0];



Beep[]
