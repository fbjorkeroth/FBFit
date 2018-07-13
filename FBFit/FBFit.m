(* ::Package:: *)

BeginPackage["FBFit`",{"FBFit`CalculateParameters`","FBFit`BestFitsAndErrors`","FBFit`SetOptions`","FBFit`Analysis`"}];

FBLoadModel::usage="Initialises the Yukawa matrices and associated variables.";
FBSetSeed::usage="Choose the seed input values for the MCMC.";
FBMonteCarlo::usage="Runs a Monte Carlo chain of length nMCMC.";

Begin["`Private`"];

(* ::Global variables:: *)

startBounds = Global`StartBounds;
isReal = Global`IsReal;
isPhase = Global`IsPhase;

(* ::Function options:: *)

Options[FBSetSeed]={"SeedSignFlip"->False,"SeedSmear"->False};
Options[FBMonteCarlo]={"Model"->"MSSM","ScaleMu"->1*^12,"TanB"->5.,"EtaB"->0.,"VaryAcceptance"->True,
	"BurnIn"->0,"SigmaGetNew"->0.01,"SaveOutput"->True,"ThinningSaveFile"->1,"Sector"->"All"};

(* ::Public functions:: *)

FBLoadModel[filename_]:=Module[{f=filename,path},
	path=FileNameJoin@{NotebookDirectory[],f};
	If[FileExistsQ[path],
		Get[path];
		Print["FBLoadModel: model loaded from ",f],
		Print["FBLoadModel: file not found!"]
	];
	Return[0]
];

FBSetSeed[seed_:Null,OptionsPattern[]]:=Module[{theta},
	If[seed === Null, SeedRandom[], SeedRandom[seed]];
	theta=N[RandomReal/@startBounds];
	If[OptionValue["SeedSignFlip"],theta=flip[theta]];
	If[OptionValue["SeedSmear"],theta=smear[theta]];
	theta
];

FBMonteCarlo[nMCMC_,theta_,OptionsPattern[]]:=Module[{t=theta,sigma,b,dbf,derr,time,l,r,tnew,lnew,alpha,meanAlpha=1.,rdata,ralpha,prog=0},
	Print["FBMonteCarlo: running fit.."];
	
	If[definedQ[OptionValue["Sector"]]==False,Print["FBMonteCarlo: global variables not defined! Quitting kernel for safety."];Quit[]];
	
	sigma=OptionValue["SigmaGetNew"];	
	b=checkBurnIn[OptionValue["BurnIn"],nMCMC];
	
	dbf=FBGetDataBestFit[];
	derr=FBGetDataErrors[];
	
	time=Now;
	l=likelihood[t,dbf,derr];
	
	Print[ProgressIndicator[Dynamic[N[prog/nMCMC]]]];
	r=Reap[Do[
		tnew=getNewTheta[t,sigma];
		lnew=likelihood[tnew,dbf,derr];
		
		alpha=Min[lnew/l,1]; (* Acceptance ratio *)
		
		If[RandomReal[]<alpha,{t,l}={tnew,lnew}]; (* Selects among old and new links *)
		
		If[OptionValue["VaryAcceptance"]==True&&(0.001<sigma<0.1),{sigma,meanAlpha}=updateSigma[sigma,meanAlpha,alpha,n]]; 
		If[n>b,
			Sow[Flatten[{t,-2Log[l]}],"Data"];
			Sow[alpha,"Alpha"]
		];
		prog++
	,{n,nMCMC}
	],{"Data","Alpha"}];

	rdata = r[[2,1,1]];
	ralpha = r[[2,2,1]];
	
	time=DateDifference[time, Now, {"Hour", "Minute"}];
	If[OptionValue["SaveOutput"],
		Export["rundata.txt",rdata[[;;;;OptionValue["ThinningSaveFile"]]],"Table"];
		Export["acceptance.log",ralpha[[;;;;OptionValue["ThinningSaveFile"]]],"List"];
		Export["time.log",ToString@time,"Text"]
	];
	Print["FBMonteCarlo: fit complete!"];
	Return[r[[2]]]
];

(* ::Internal functions:: *)

definedQ[sec_]:=Module[{q,l,inp,is,checklist},
	q={ValueQ[Global`Yu],ValueQ[Global`Yd]};
	l={ValueQ[Global`Mnu],ValueQ[Global`Ye]};
	inp={ValueQ[Global`InputVariables],ValueQ[Global`StartBounds],ValueQ[Global`IsReal],ValueQ[Global`IsPhase]};
	is={ValueQ[Global`IsQuark],ValueQ[Global`IsLepton]};
	checklist=Switch[sec,
		"Q",Join[q,inp],
		"L",Join[l,inp],
		_,Join[q,l,inp,is]
	];
	And@@checklist
];

likelihood[theta_,databestfit_,dataerrors_]:=Module[{ca,pu,chi2},
	ca=FBGetPhysicalParameters[theta];
	pu=FBGetPulls[ca,databestfit,dataerrors];
	chi2=FBChiSq[pu];
	Exp[-chi2/2]
];

checkBurnIn[burnin_,n_]:=If[burnin>n,
	Print["FBMonteCarlo: burn-in larger than chain length! Setting to zero.."];0,
	burnin
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
