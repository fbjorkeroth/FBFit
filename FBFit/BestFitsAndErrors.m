(* ::Package:: *)

BeginPackage["FBFit`BestFitsAndErrors`"];

FBLoadBestFitsAndErrors::usage="Generates masses and mixing parameters (in MSSM, defines functions in terms of tan \[Beta] and etaB for a given value of MSUSY).";
FBGetDataBestFit::usage="Generates list of experimental best fit points.";
FBGetDataErrors::usage="Generates list of experimental errors.";

Begin["`Private`"];

(* ::Function options:: *)

Options[FBLoadBestFitsAndErrors]={"Model"->"SM","MSUSY"->1,"ScaleMu"->100,"NeutrinoOrdering"->"Normal","UniversalError"->Null};
Options[FBGetDataBestFit]={"Model"->"SM","ScaleMu"->100,"TanB"->5.,"EtaB"->0.,"Sector"->"All"};
Options[FBGetDataErrors]={"Model"->"SM","ScaleMu"->100,"TanB"->5.,"EtaB"->0.,"Sector"->"All","ExcludeParameters"->{}};

(* ::Public functions:: *)

FBLoadBestFitsAndErrors[opts:OptionsPattern[]]:=Module[{model},
	model=OptionValue["Model"];
	Print["FBLoadBestFitsAndErrors: extracting Yukawa couplings and mixing parameters..."];
	Switch[model,
	"MSSM",loadBestFitsAndErrorsMSSM[OptionValue["MSUSY"],OptionValue["NeutrinoOrdering"]],
	"SM",loadBestFitsAndErrorsSM[OptionValue["ScaleMu"],OptionValue["NeutrinoOrdering"],OptionValue["UniversalError"]],
	_,printBadModel[]
	]
];

FBGetDataBestFit[OptionsPattern[]]:=Module[{model,sec},
	model=OptionValue["Model"];
	sec=OptionValue["Sector"];
	Switch[model,
	"MSSM",getDataBestFitMSSM[sec,OptionValue["TanB"],OptionValue["EtaB"]],
	"SM",getDataBestFitSM[sec,OptionValue["ScaleMu"]],
	_,printBadModel[]
	]
];

FBGetDataErrors[OptionsPattern[]]:=Module[{model,sec,err},
	model=OptionValue["Model"];
	sec=OptionValue["Sector"];
	err=Switch[model,
	"MSSM",getDataErrorsMSSM[sec,OptionValue["TanB"],OptionValue["EtaB"]],
	"SM",getDataErrorsSM[sec,OptionValue["ScaleMu"]],
	_,printBadModel[]
	];	
	ReplacePart[err,#->10^8.&/@OptionValue["ExcludeParameters"]]
];

(* ::Internal functions:: *)

printBadModel[]:=Module[{},Print["Model not supported. Quitting kernel for safety.."];Quit[]];

loadBestFitsAndErrorsSM[mu_,ordering_,univerror_]:=Module[{vHiggs=174},
	{theta12q,theta13q,theta23q,deltaq}={13.026,3.8*^-3/Degree,4.4*^-2/Degree,69.215};(* Cabibbo angle and CP phase *)
	{errtheta12q,errtheta13q,errtheta23q,errdeltaq}={0.041,0.036 theta13q,0.016 theta23q,3.095};

	{yu,yc,yt}={0.58*^-3,0.281,86.7}/vHiggs;
	{erryu,erryc,erryt}={0.24*^-3,0.042,4.0}/vHiggs;

	{yd,ys,yb}={1.34*^-3,26.*^-3,1.21}/vHiggs;
	{erryd,errys,erryb}={6.0*^-4,8*^-3,0.05}/vHiggs;

	{ye,ymu,ytau}={4.90856087*^-4,0.103622931,1.76167}/vHiggs;
	{errye,errymu,errytau}=0.006{ye,ymu,ytau};

(*	(* Lepton data from NuFit 3.2 *)
	{theta12l,theta13l,theta23l,deltal}={33.62,8.54,47.2,234};
	{errtheta12l,errtheta13l,errtheta23l,errdeltal}={0.78,0.15,3.9,43};
	
	{dm21,dm31}={7.40*^-5,2.494*^-3};
	{errdm21,errdm31}={0.21*^-5,0.033*^-3};*)	
	
	If[mu=="MZ",
		{theta12q,theta13q,theta23q,deltaq}={0.22735,3.64*^-3,4.208*^-2,1.208}/Degree;(* Cabibbo angle and CP phase *)
		{errtheta12q,errtheta13q,errtheta23q,errdeltaq}={0.00072,0.13*^-3,0.064*^-2,0.054}/Degree;

		{yu,yc,yt}={7.4*^-6,3.60*^-3,0.9861};
		{erryu,erryc,erryt}={3.0*^-6,0.11*^-3,0.0086};

		{yd,ys,yb}={1.58*^-5,3.12*^-4,1.639*^-2};
		{erryd,errys,erryb}={0.23*^-5,0.17*^-4,0.015*^-2};

		{ye,ymu,ytau}={4.90856087*^-4,0.103622931,1.76167}/vHiggs;
		{errye,errymu,errytau}=0.006{ye,ymu,ytau};
	];
	
	(* Data from NuFit 4.0 *)
	Switch[ordering,
		"Normal",
			{theta12l,theta13l,theta23l,deltal}={33.82,8.61,49.7,217};
			{errtheta12l,errtheta13l,errtheta23l,errdeltal}={0.78,0.13,1.1,40};

			{dm21,dm31}={7.39*^-5,2.525*^-3};
			{errdm21,errdm31}={0.21*^-5,0.033*^-3};
		,
		"Inverted",
			{theta12l,theta13l,theta23l,deltal}={33.82,8.65,49.7,280};
			{errtheta12l,errtheta13l,errtheta23l,errdeltal}={0.78,0.13,1.0,28};
							
			{dm21,dm31}={7.39*^-5,-2.512*^-3};
			{errdm21,errdm31}={0.21*^-5,0.034*^-3};
	];
	
	If[univerror!=Null,
		{errtheta12q,errtheta13q,errtheta23q,errdeltaq}=univerror{theta12q,theta13q,theta23q,deltaq};
		{erryu,erryc,erryt,erryd,errys,erryb}=univerror{yu,yc,yt,yd,ys,yb};
		{errye,errymu,errytau,errdm21,errdm31}=univerror{ye,ymu,ytau,dm21,dm31}
	];
	
	Return[0];
];

loadBestFitsAndErrorsMSSM[MSUSY_,ordering_]:=Module[{data(*,vHiggs=174*)},
	SetOptions[Interpolation,InterpolationOrder->1];

	data=importDataFile[#,MSUSY]&/@{
		"reduced_theta13","reduced_theta23",
		"reduced_yu","reduced_yd","reduced_ye",
		"reduced_yc","reduced_ys","reduced_ymu",
		"reduced_yt","reduced_yb","reduced_ytau",
		"sigma_yb","sigma_yt","sigma_ytau"
	};

	theta13q=Function[{tanBeta,etaB},10^-3. (1+etaB)Quiet@Interpolation[data[[1]]][tanBeta,etaB]/Degree];
	theta23q=Function[{tanBeta,etaB},10^-2. (1+etaB)Quiet@Interpolation[data[[2]]][tanBeta,etaB]/Degree];
	yu=Function[{tanBeta,etaB},10^-6 Csc[ArcTan[tanBeta]]Quiet@Interpolation[data[[3]]][tanBeta,etaB]];
	yd=Function[{tanBeta,etaB},10^-5 Sec[ArcTan[tanBeta]]Quiet@Interpolation[data[[4]]][tanBeta,etaB]];
	ye=Function[{tanBeta,etaB},10^-6 Sec[ArcTan[tanBeta]]Quiet@Interpolation[data[[5]]][tanBeta,etaB]];
	yc=Function[{tanBeta,etaB},10^-3 Csc[ArcTan[tanBeta]]Quiet@Interpolation[data[[6]]][tanBeta,etaB]];
	ys=Function[{tanBeta,etaB},10^-4 Sec[ArcTan[tanBeta]]Quiet@Interpolation[data[[7]]][tanBeta,etaB]];
	ymu=Function[{tanBeta,etaB},10^-4 Sec[ArcTan[tanBeta]]Quiet@Interpolation[data[[8]]][tanBeta,etaB]];
	yt=Function[{tanBeta,etaB},Csc[ArcTan[tanBeta]]Quiet@Interpolation[data[[9]]][tanBeta,etaB]];
	yb=Function[{tanBeta,etaB},10^-2 (1+etaB)^-1 Sec[ArcTan[tanBeta]]Quiet@Interpolation[data[[10]]][tanBeta,etaB]];
	ytau=Function[{tanBeta,etaB},10^-2 Sec[ArcTan[tanBeta]]Quiet@Interpolation[data[[11]]][tanBeta,etaB]];
	
	sigmayt=Function[{tanBeta,etaB},Quiet@Interpolation[data[[12]]][tanBeta,etaB]];
	sigmayb=Function[{tanBeta,etaB},Quiet@Interpolation[data[[13]]][tanBeta,etaB]];
	sigmaytau=Function[{tanBeta,etaB},Quiet@Interpolation[data[[14]]][tanBeta,etaB]];
	
	errtheta13q=Function[{tanBeta,etaB},0.036theta13q[tanBeta,etaB]];
	errtheta23q=Function[{tanBeta,etaB},0.016theta23q[tanBeta,etaB]];
	erryu=Function[{tanBeta,etaB},0.31yu[tanBeta,etaB]];
	erryd=Function[{tanBeta,etaB},0.11yd[tanBeta,etaB]];
	errye=Function[{tanBeta,etaB},0.006ye[tanBeta,etaB]];
	erryc=Function[{tanBeta,etaB},0.035yc[tanBeta,etaB]];
	errys=Function[{tanBeta,etaB},0.054ys[tanBeta,etaB]];
	errymu=Function[{tanBeta,etaB},0.006ymu[tanBeta,etaB]];
	erryt=Function[{tanBeta,etaB},sigmayt[tanBeta,etaB]yt[tanBeta,etaB]];
	erryb=Function[{tanBeta,etaB},sigmayb[tanBeta,etaB]yb[tanBeta,etaB]];
	errytau=Function[{tanBeta,etaB},sigmaytau[tanBeta,etaB]ytau[tanBeta,etaB]];
	
	(*{theta12q,deltaq}={13.026,69.215};
	{theta12l,theta13l,theta23l}={33.57,8.46,41.75};
	{dm21,dm31}={7.51*^-5,2.524*^-3};
	deltal=257;

	{errtheta12q,errdeltaq}={0.041,3.095};
	{errtheta12l,errtheta13l,errtheta23l}={0.76,0.15,1.35};
	{errdm21,errdm31}={0.18*^-5,0.040*^-3};
	errdeltal=55;*)
	
	(* Data from NuFit 4.0 *)
	Switch[ordering,
		"Normal",
			{theta12l,theta13l,theta23l,deltal}={33.82,8.61,49.7,217};
			{errtheta12l,errtheta13l,errtheta23l,errdeltal}={0.78,0.13,1.1,40};

			{dm21,dm31}={7.39*^-5,2.525*^-3};
			{errdm21,errdm31}={0.21*^-5,0.033*^-3};
		,
		"Inverted",
			{theta12l,theta13l,theta23l,deltal}={33.82,8.65,49.7,280};
			{errtheta12l,errtheta13l,errtheta23l,errdeltal}={0.78,0.13,1.0,28};
							
			{dm21,dm31}={7.39*^-5,-2.512*^-3};
			{errdm21,errdm31}={0.21*^-5,0.034*^-3};
	];
];

thisPath=DirectoryName[$InputFileName];
importDataFile[yuk_,MSUSY_]:=Module[{dDir,endName,fileName},
	endName="MSUSY="<>ToString[MSUSY]<>"TeV/";
	fileName=yuk<>".dat";
	dDir=FileNameJoin@{
		thisPath,
		"MixingParameterTools",
		"RunningParameters",
		"data",
		endName,
		fileName
	};
	Import[dDir,"Table"]
];



getDataBestFitSM[sec_,mu_]:=getDataBestFitSM[sec,mu]=Module[{q,l},
	q={theta12q,theta13q,theta23q,deltaq,yu,yc,yt,yd,ys,yb};
	l={theta12l,theta13l,theta23l,deltal,dm21,dm31,ye,ymu,ytau};
	Switch[sec,"Q",q,"L",l,_,Join[q,l]]
];

getDataErrorsSM[sec_,mu_]:=getDataErrorsSM[sec,mu]=Module[{q,l},
	q={errtheta12q,errtheta13q,errtheta23q,errdeltaq,erryu,erryc,erryt,erryd,errys,erryb};
	l={errtheta12l,errtheta13l,errtheta23l,errdeltal,errdm21,errdm31,errye,errymu,errytau};
	Switch[sec,"Q",q,"L",l,_,Join[q,l]]
];

getDataBestFitMSSM[sec_,tanBeta_,etaB_]:=getDataBestFitMSSM[sec,tanBeta,etaB]=Module[{q,l},
	q={theta12q,theta13q[tanBeta,etaB],theta23q[tanBeta,etaB],deltaq,
		yu[tanBeta,etaB],yc[tanBeta,etaB],yt[tanBeta,etaB],
		yd[tanBeta,etaB],ys[tanBeta,etaB],yb[tanBeta,etaB]};
	l={theta12l,theta13l,theta23l,deltal,
		dm21,dm31,
		ye[tanBeta,etaB],ymu[tanBeta,etaB],ytau[tanBeta,etaB]};
	Switch[sec,"Q",q,"L",l,_,Join[q,l]]
];

getDataErrorsMSSM[sec_,tanBeta_,etaB_]:=getDataErrorsMSSM[sec,tanBeta,etaB]=Module[{q,l},
	q={errtheta12q,errtheta13q[tanBeta,etaB],errtheta23q[tanBeta,etaB],errdeltaq,
		erryu[tanBeta,etaB],erryc[tanBeta,etaB],erryt[tanBeta,etaB],
		erryd[tanBeta,etaB],errys[tanBeta,etaB],erryb[tanBeta,etaB]};
	l={errtheta12l,errtheta13l,errtheta23l,errdeltal,
		errdm21,errdm31,
		errye[tanBeta,etaB],errymu[tanBeta,etaB],errytau[tanBeta,etaB]};
	Switch[sec,"Q",q,"L",l,_,Join[q,l]]
];

End[];
EndPackage[];
