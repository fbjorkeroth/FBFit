(* ::Package:: *)

BeginPackage["FBFit`SetOptions`"];

FBSetOptions::usage="Sets all function options to values defined by the user.";

Begin["Private`"];

(* ::Function options:: *)

Options[FBSetOptions]={
	"SaveOutput"->True,
	"Analysis"->False,
	"Model"->"SM",
	"MSUSY"->1,
	"ScaleMu"->100,
	"TanB"->5.,
	"EtaB"->0.,
	"BurnIn"->0,
	"VarySigma"->False,
	"ExcludeParameters"->{},
	"SeedSignFlip"->False,
	"SeedSmear"->False,
	"SigmaGetNew"->0.01,
	"ThinningSaveFile"->1,
	"Sector"->"All",
	"MinAcceptance"->0,
	"NeutrinoOrdering"->"Normal",
	"UniversalError"->Null
};

(* ::Public functions:: *)

FBSetOptions[opts:OptionsPattern[]]:=Module[{save,an,m,ms,scmu,tb,eb,bu,va,ep,ssf,ss,sigma,tsf,sec,mac,ord,ue},
	{save,an,m,ms,scmu,tb,eb,bu,va,ep,ssf,ss,sigma,tsf,sec,mac,ord,ue}=OptionValue[#]&/@{
		"SaveOutput",
		"Analysis",
		"Model",
		"MSUSY",
		"ScaleMu",
		"TanB",
		"EtaB",
		"BurnIn",
		"VarySigma",
		"ExcludeParameters",
		"SeedSignFlip",
		"SeedSmear",
		"SigmaGetNew",
		"ThinningSaveFile",
		"Sector",
		"MinAcceptance",
		"NeutrinoOrdering",
		"UniversalError"
		};
		
	SetOptions[FBFit`FBSetSeed,{"SeedSignFlip"->ssf,"SeedSmear"->ss}];
	SetOptions[FBFit`FBMonteCarlo,{"Model"->m,"ScaleMu"->scmu,"TanB"->tb,"EtaB"->eb,"BurnIn"->bu,"VarySigma"->va,
		"SigmaGetNew"->sigma,"SaveOutput"->save,"ThinningSaveFile"->tsf,"Sector"->sec,"MinAcceptance"->mac}];

	SetOptions[FBFit`CalculateParameters`FBGetPhysicalParameters,{"Sector"->sec}];

	SetOptions[FBFit`BestFitsAndErrors`FBLoadBestFitsAndErrors,{"Model"->m,"MSUSY"->ms,"ScaleMu"->scmu,"NeutrinoOrdering"->ord,"UniversalError"->ue}];
	SetOptions[FBFit`BestFitsAndErrors`FBGetDataBestFit,{"Model"->m,"ScaleMu"->scmu,"TanB"->tb,"EtaB"->eb,"Sector"->sec}];
	SetOptions[FBFit`BestFitsAndErrors`FBGetDataErrors,{"Model"->m,"ScaleMu"->scmu,"TanB"->tb,"EtaB"->eb,
		"Sector"->sec,"ExcludeParameters"->ep}];

	SetOptions[FBFit`Analysis`FBPrintInput,{"Model"->m,"TanB"->tb,"EtaB"->eb}];
	SetOptions[FBFit`Analysis`FBPrintOutput,{"Model"->m,"TanB"->tb,"EtaB"->eb,"Sector"->sec}];
	SetOptions[FBFit`Analysis`FBPlotPulls,{"Model"->m,"TanB"->tb,"EtaB"->eb,"Sector"->sec}];

	SetOptions[FBSetOptions,opts];
	Print["FBSetOptions: fit specs set: ",ToString[{opts}]];
	Which[
		an==True,
			Print["FBSetOptions: options set for analysis."],
		save==True,
			Print["FBSetOptions: building directories..."];
			buildDataDirectories[];
			SetOptions[FBSetOptions,"Analysis"->True];
			Export["variables.mx",Options[FBSetOptions]],
		save==False,
			Print["FBSetOptions: no data will be exported!"]
	];
	{opts}
];


(* ::Internal functions:: *)

buildDataDirectories[]:=Module[{dataHomeDir,runName,runComments,runDir,validDir=False,ack},
	dataHomeDir=FileNameJoin@{NotebookDirectory[],"data"};
	If[DirectoryQ[dataHomeDir]==False,
		CreateDirectory[dataHomeDir];
		CreateDialog[{TextCell["Created new directory for storing data."],DefaultButton[DialogReturn[ack=True]]},Modal->True],
		ack=True (* ack: requires user to acknowledge that a data directory has been created *)
	];
	
	WaitUntil[ack==True];
	{runName,runComments}=inputName[];
	While[validDir==False,
		runDir=FileNameJoin@{dataHomeDir,runName};
		Which[
		runName!="" && DirectoryQ[runDir]==False,
			validDir=True,
		runName=="",
			{runName,runComments}=inputName["blank"];
			validDir=False,
		DirectoryQ[runDir],
			{runName,runComments}=inputName["overwrite"];
			validDir=False
		]
	];

	CreateDirectory[runDir];
	SetDirectory[runDir];
	Print[
		"Working directory set to ",runDir,
		"\nAll data and parameters will be saved here."];
	Export["comments.txt",runComments];
	Return[runDir]
];


inputName[invalidFlag_:""]:=DialogInput[{name="",comments=""},
	Column[
	{
	Which[
		invalidFlag=="blank","Run name can't be null! Try again: ",
		invalidFlag=="overwrite","Directory already exists! Try again: ",
		invalidFlag=="","Directory name for this run: "
	],
	InputField[Dynamic[name],String],
	"Comments (optional): ",
	InputField[Dynamic[comments],String],
	Button["OK",DialogReturn[{name,comments}],ImageSize->Automatic]
	}
	]
];


End[];
EndPackage[];
