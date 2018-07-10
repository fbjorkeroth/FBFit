(* ::Package:: *)

(* ::Title:: *)
(*Simple MCMC*)


(* ::Section:: *)
(*Setup*)


AppendTo[$Path,ParentDirectory[NotebookDirectory[]]];
Get["FBFit`"];

FBLoadModel["models/model.m"];

nMCMC=100;(* Number of generated sets. Roughly 1 min per 40k on iMac (4 x i7) *)
nThin=Max[Round[nMCMC/50000],1];(* Take 1/nThin of generated points *)

FBSetOptions[
	"Model"->"SM",
	"BurnIn"->0,
	"SaveOutput"->False,
	"ExcludeParameters"->Range[11,19],
	"Thinning"->1
];
FBLoadBestFitsAndErrors[];


(* ::Section:: *)
(*Run fit*)


Now
\[Theta]0=FBSetSeed[];
l=FBMonteCarlo[nMCMC,\[Theta]0];
Now


(* ::Section:: *)
(*Endcap*)


Beep[]
