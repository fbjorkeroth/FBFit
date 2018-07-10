(* ::Package:: *)

BeginPackage["FBFit`SetOptions`"];

FBSetOptions::usage="Sets all function options to values defined by the user.";

Begin["Private`"];

(* ::Function options:: *)

Options[FBSetOptions]={
	"SaveOutput"->True,
	"Analysis"->False,
	"Model"->"MSSM",
	"MSUSY"->1,
	"ScaleMu"->1*^12,
	"TanB"->5.,
	"EtaB"->0.,
	"BurnIn"->0,
	"VaryAcceptance"->True,
	"ExcludeParameters"->{},
	"SeedSignFlip"->False,
	"SeedSmear"->False,
	"SigmaGetNew"->0.01,
	"ThinningSaveFile"->1,
	"Thinning"->1
};

(* ::Public functions:: *)

FBSetOptions[opts:OptionsPattern[]]:=Module[{save,an,m,ms,scmu,tb,eb,bu,va,ep,ssf,ss,sigma,tsf},
	{save,an,m,ms,scmu,tb,eb,bu,va,ep,ssf,ss,sigma,tsf}=OptionValue[#]&/@{
		"SaveOutput",
		"Analysis",
		"Model",
		"MSUSY",
		"ScaleMu",
		"TanB",
		"EtaB",
		"BurnIn",
		"VaryAcceptance",
		"ExcludeParameters",
		"SeedSignFlip",
		"SeedSmear",
		"SigmaGetNew",
		"ThinningSaveFile"
		};
		
	SetOptions[FBFit`FBSetSeed,{"SeedSignFlip"->ssf,"SeedSmear"->ss}];
	SetOptions[FBFit`FBMonteCarlo,{"Model"->m,"ScaleMu"->scmu,"TanB"->tb,"EtaB"->eb,"BurnIn"->bu,
			"VaryAcceptance"->va,"ExcludeParameters"->ep,"SigmaGetNew"->sigma,"SaveOutput"->save,"ThinningSaveFile"->tsf}];

	SetOptions[FBFit`BestFitsAndErrors`FBLoadBestFitsAndErrors,{"Model"->m,"MSUSY"->ms,"ScaleMu"->scmu}];
	SetOptions[FBFit`BestFitsAndErrors`FBGetDataBestFit,{"Model"->m,"ScaleMu"->scmu,"TanB"->tb,"EtaB"->eb}];
	SetOptions[FBFit`BestFitsAndErrors`FBGetDataErrors,{"Model"->m,"ScaleMu"->scmu,"TanB"->tb,"EtaB"->eb}];

	SetOptions[FBFit`Analysis`FBPrintInput,{"Model"->m,"TanB"->tb,"EtaB"->eb,"ExcludeParameters"->ep}];
	SetOptions[FBFit`Analysis`FBPrintOutput,{"Model"->m,"TanB"->tb,"EtaB"->eb}];
	SetOptions[FBFit`Analysis`FBPlotPulls,{"Model"->m,"TanB"->tb,"EtaB"->eb,"ExcludeParameters"->ep}];

	SetOptions[FBSetOptions,opts];
	Print["FBSetOptions: model specs set: ",ToString[{opts}]];
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
