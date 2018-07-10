(* ::Package:: *)

BeginPackage["FBFit`",{"FBFit`CalculateParameters`","FBFit`BestFitsAndErrors`","FBFit`SetOptions`","FBFit`Analysis`"}];

FBLoadModel::usage="Initialises the Yukawa matrices and associated variables.";
FBSetSeed::usage="Choose the seed input values for the MCMC.";
FBMonteCarlo::usage="Runs a Monte Carlo chain of length nMCMC.";

Begin["`Private`"];

(* ::Function options:: *)

Options[FBSetSeed]={"SeedSignFlip"->False,"SeedSmear"->False};
Options[FBMonteCarlo]={"Model"->"MSSM","ScaleMu"->1*^12,"TanB"->5.,"EtaB"->0.,"ExcludeParameters"->{},
	"VaryAcceptance"->True,"BurnIn"->0,"SigmaGetNew"->0.01,"SaveOutput"->True,"Thinning"->1};

(* ::Global variables:: *)

Yu=Global`Yu;
Yd=Global`Yd;
mnu=Global`mnu;
Ye=Global`Ye;

inputVariables = Global`inputVariables;
inLabels = Global`inLabels;
startBounds = Global`startBounds;

isReal = Global`isReal;
isPhase = Global`isPhase;
isQuark = Global`isQuark;
isLepton = Global`isLepton;

(* ::Public functions:: *)

FBLoadModel[filename_]:=Module[{f=filename},
	Get[FileNameJoin@{NotebookDirectory[],f}];
	Print["FBLoadModel: model loaded from ",f]
];


FBSetSeed[OptionsPattern[]]:=Module[{theta},
	theta=N[RandomReal/@startBounds];
	If[OptionValue["SeedSignFlip"],theta=flip[theta]];
	If[OptionValue["SeedSmear"],theta=smear[theta]];
	theta
];

FBMonteCarlo[nMCMC0_,theta0_,OptionsPattern[]]:=Module[{nMCMC=nMCMC0,t=theta0,sigma,dbf,derr,time,l,r,tnew,lnew,alpha,meanAlpha=1.,rdata,ralpha},
	
	sigma=OptionValue["SigmaGetNew"];
	
	dbf=FBGetDataBestFit[];
	derr=ReplacePart[FBGetDataErrors[],#->10^8.&/@OptionValue["ExcludeParameters"]];
	
	time=Now;
	l=likelihood[t,dbf,derr];
	r=Reap[Do[(* Calculates the new link in the MCMC chain. *)
		tnew=getNewTheta[t,sigma]; (* Get input parameters for new link *)
		lnew=likelihood[tnew,dbf,derr];
		
		alpha=Min[lnew/l,1]; (* Acceptance ratio *)
		If[RandomReal[]<alpha,{t,l}={tnew,lnew}]; (* Selects among old and new links *)

		If[OptionValue["VaryAcceptance"]==True,{sigma,meanAlpha}=updateSigma[sigma,meanAlpha,alpha,n]]; 
		(* Functionality for progressively updating the Gaussian width. *)

		If[n>OptionValue["BurnIn"],
			Sow[Flatten[{t,-2Log[l]}],"Data"];
			Sow[alpha,"Alpha"]
		](* Adds the link to the output chain *)
	,{n,nMCMC}
	],{"Data","Alpha"}];

	rdata = r[[2,1,1]];
	ralpha = r[[2,2,1]];
	
	time=DateDifference[time, Now, {"Hour", "Minute"}];
	If[OptionValue["SaveOutput"],
		Export["rundata.txt",rdata[[;;;;OptionValue["Thinning"]]],"Table"];
		Export["acceptance.log",ralpha[[;;;;OptionValue["Thinning"]]],"List"];
		Export["time.log",ToString@time,"Text"]
	];
	Return[r]
];

(* ::Internal functions:: *)

likelihood[theta_,databestfit_,dataerrors_]:=Module[{m,ca,pu,chi2},
	m={Yu,Yd,mnu,Ye}/.Thread[inputVariables->theta];
	ca=FBCalculateParameters@@m;
	pu=FBGetPulls[ca,databestfit,dataerrors];
	chi2=FBChiSq[pu];
	Exp[-chi2/2]
];

printBadModel[]:=Module[{},Print["Model not supported. Quitting kernel for safety.."];Quit[]];

flip[theta_]:=Table[RandomChoice@{1,-1}theta[[i]],{i,Length[theta]}];
smear[theta_,f_:0.2]:=RandomReal[{(1-f)#,(1+f)#}]&/@theta;
(* Functions allowing you to randomly change the signs of or smear the starting value of theta *)

getNewTheta[theta_,sigma_,fixlist_:{}]:=Module[{tnew,select},
	select[list_]:=DeleteCases[list,n_/;MemberQ[fixlist,n]]; (* Fixes parameters given by fixlist *)
	tnew=MapAt[findYukawa[#,sigma]&,theta,{#}&/@select[isReal]]; (* Updates real parameters *)
	MapAt[findPhase[#,sigma]&,tnew,{#}&/@select[isPhase]] (* Updates phase parameters *)
];

findYukawa[t_,sigma_,floor_:10^-7]:=RandomVariate[NormalDistribution[t,Max[Abs[sigma t],floor]]];
findPhase[t_,sigma_]:=Mod[RandomVariate[NormalDistribution[t,sigma]],2\[Pi]];

updateSigma[sigma_,meanalpha_,alpha_,n_]:=Module[{newmean},
	newmean=(n-1)meanalpha/n+alpha/n;
	Return[{
	Which[
		newmean<0.3,0.99sigma,
		newmean>0.5,1.01sigma,
		0.3<=newmean<=0.5,sigma
	],
	newmean
	}]
];

End[];
EndPackage[];
