(* ::Package:: *)

(* ::Title:: *)
(*Simple MCMC*)


AppendTo[$Path,ParentDirectory[NotebookDirectory[]]];
Get["FBFit`"];

FBLoadModel["models/model.m"];

nMCMC=1000;

FBSetOptions[
	"Model"->"SM",
	"SaveOutput"->False
];

FBLoadBestFitsAndErrors[];

\[Theta]0=FBSetSeed[];

l=FBMonteCarlo[nMCMC,\[Theta]0];

Print["Done!"];



