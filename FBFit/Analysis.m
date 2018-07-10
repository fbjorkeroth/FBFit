(* ::Package:: *)

BeginPackage["FBFit`Analysis`",{"FBFit`CalculateParameters`","FBFit`BestFitsAndErrors`"}];
Quiet[Needs["Combinatorica`"],General::compat];

FBImportFrom::usage="Sets the data directory to import from.";
FBExtractBestInput::usage="Extracts from the dataset the input parameters for the best fit point.";

FBPrintInput::usage="prints the best fit \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\) and input parameters.";
FBPrintOutput::usage="prints the best fit physical parameters.";
FBPlotPulls::usage="plots a bar chart of the pulls of each output parameter.";

FBPhysicalParameters::usage="Outputs the physical parameters (Yukawa eigenvalues, mixing parameters) for a given set of input data. Option exists to thin this data for quicker evaluation.";
FBCredibleInterval::usage="Calculates the N% credible interval for physical (Yukawa, mixing) parameters for a given input data set.";
FBPlotHistogram::usage="Plots histogram of physical output data.";

FBChopDataFraction::usage="Chops off the first Nth fraction of rows in a data set, where 0<N<1.";
FBAverageRejection::usage="Calculates the average acceptance rate by counting successive elements of output table with the same \[Chi]2. Essentially defunct now.";


Begin["Private`"];

(* ::Function options:: *)

Options[FBPrintInput]={"Model"->"MSSM","TanB"->5.,"EtaB"->0.,"ExcludeParameters"->{}};
Options[FBPrintOutput]={"Model"->"MSSM","TanB"->5.,"EtaB"->0.};
Options[FBPlotPulls]={"Model"->"MSSM","TanB"->5.,"EtaB"->0.,"ExcludeParameters"->{}};

Options[FBPhysicalParameters]={"Thinning"->1};
Options[FBCredibleInterval]={"Thinning"->1,"SigmaSpan"->3,"PixelDensity"->100,"Plot"->True};
Options[FBPlotHistogram]={"Thinning"->1,"ImageSize"->300,"Bins"->30};

(* ::Global variables:: *)

Yu = Global`Yu;
Yd = Global`Yd;
mnu = Global`mnu;
Ye = Global`Ye;

inputVariables = Global`inputVariables;
inLabels = Global`inLabels;
startBounds = Global`startBounds;

isReal = Global`isReal;
isPhase = Global`isPhase;
isQuark = Global`isQuark;
isLepton = Global`isLepton;

(* ::Public functions:: *)

FBImportFrom[runName_]:=Module[{rn=runName,dir},
	dir=FileNameJoin@{NotebookDirectory[],"data",rn};
	If[DirectoryQ[dir]==False,
		Print["FBImportFrom: data directory ",dir," not found! Please select another run. Quitting kernel for safety.."];
		Quit[]
	];
	SetDirectory[dir];
	Print["FBImportFrom: directory set to ",dir]
];

FBExtractBestInput[inputdata_]:=Module[{data=inputdata,lRow,chisq,best},
	lRow=Length[data[[1]]];
	chisq=data[[;;,lRow]];
	best=Position[chisq,Min[chisq]][[1,1]];
	data[[best,;;-2]]
];

FBPrintInput[theta_,OptionsPattern[]]:=Module[{t=theta,databestfit,dataerrors,titles,chisq,table},
	databestfit=FBGetDataBestFit[];
	dataerrors=ReplacePart[FBGetDataErrors[],#->10^8.&/@OptionValue["ExcludeParameters"]];
	titles={"Parameter","Value"};
	chisq=FBChiSq[FBGetPulls[FBCalculateParameters[Yu,Yd,mnu,Ye],databestfit,dataerrors]]/.Thread[inputVariables->t];
	table=MatrixForm@Prepend[Transpose@{inLabels[[#]],t[[#]]},titles]&/@{isQuark,isLepton};
	Print[
		"\!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\): ",chisq,
		"\nInput: (quarks) ",table[[1]],"\t(leptons) ",table[[2]]
	]
];

FBPrintOutput[theta_,OptionsPattern[]]:=Module[{t=theta,databestfit,calc,titles,table,x},
	databestfit=FBGetDataBestFit[];
	calc=FBCalculateParameters[Yu,Yd,mnu,Ye]/.Thread[inputVariables->t];
	titles={"Parameter","Data","Model"};
	table=MatrixForm@Prepend[Transpose@{outLabels,databestfit,calc},titles];
	x=ScientificForm[table,4,ExponentFunction->(If[-2<#<3,Null,#]&)];
	Print["Output: ",x]
];

FBPlotPulls[theta_,OptionsPattern[]]:=Module[{t=theta,databestfit,dataerrors,calc,pulls},
	databestfit=FBGetDataBestFit[];
	dataerrors=ReplacePart[FBGetDataErrors[],#->10^8.&/@OptionValue["ExcludeParameters"]];
	calc=FBCalculateParameters[Yu,Yd,mnu,Ye]/.Thread[inputVariables->t];
	pulls=FBGetPulls[calc,databestfit,dataerrors];
	Print@BarChart[pulls,ChartLabels->outLabels,AxesLabel->"Pull",ImageSize->Large,AspectRatio->1/2,BaseStyle->FontSize->12]
];

FBPhysicalParameters[inputdata_,variable_,OptionsPattern[]]:=Module[{data=inputdata,dataThinned,n=variable},
	dataThinned=Take[data,{1,Length[data],OptionValue["Thinning"]}];
	(FBCalculateParameters[Yu,Yd,mnu,Ye]/.Thread[inputVariables->dataThinned[[#,;;-2]]])[[n]]&/@Range[Length[dataThinned]]
];


FBCredibleInterval[inputdata_,variable_,CIlevel_,OptionsPattern[]]:=Module[{data=inputdata,n=variable,level=CIlevel,l,mu,std,pdf,spacing,grid,t,plot,ci},
	l=FBPhysicalParameters[data,n,"Thinning"->OptionValue["Thinning"]];
	{mu,std,pdf}=#[l]&/@{Mean,StandardDeviation,SmoothKernelDistribution};
	spacing=10std/OptionValue["PixelDensity"];
	grid=Table[PDF[pdf,{mu-OptionValue["SigmaSpan"] std+spacing*i}],{i,OptionValue["PixelDensity"]}];
	t=findCredibilityLevel[grid,level];

	plot=Plot[PDF[pdf,x]-t,{x,mu-std*OptionValue["SigmaSpan"],mu+std*OptionValue["SigmaSpan"]},
		Mesh->{{0}},
		MeshFunctions->{#2&},
		MeshStyle->PointSize[Medium],
		PlotRange->All,
		ImageSize->Medium,
		Axes->{True,False}
	];
	ci=Sort@Cases[Normal@plot,Point[{x_,y_}]->x,Infinity];

	If[OptionValue["Plot"]==True,Print[plot]];
	Print[StringJoin["n = ",ToString[n],"; ",ToString[100level//Round],"% CI: ",ToString[ci]]];
	Return[]
];

FBPlotHistogram[inputdata_,variable_,OptionsPattern[]]:=Module[{data=inputdata,n=variable,bf,err,p,bins,h,maxH,line},
	bf=FBGetDataBestFit[][[n]];
	err=FBGetDataErrors[][[n]];
	p=FBPhysicalParameters[data,n,"Thinning"->OptionValue["Thinning"]];
	{bins, h} = HistogramList[p,OptionValue["Bins"]];
	maxH=Max[h];
	line[x_]:=Line[{{x,0},{x,maxH+2}}];
	Histogram[
		p,{bins},h&,
		Epilog->{line[bf],{Dashed,line[bf+err]},{Dashed,line[bf-err]}}
	]
];

FBChopDataFraction[inputdata_,bottomfraction_]:=Module[{data=inputdata,b=bottomfraction,n},
	n=Round[(1-b)Length[data]];
	Take[data,-n]
];


FBAverageRejection[table_,n_]:=Module[{l},
	l=Exp[-table[[;;,Length@table[[1]]]]/2];
	Count[Round@Ratios[l],1]/n//N
];

(* ::Internal functions:: *)

printBadModel[]:=Module[{},Print["Model not supported. Quitting kernel for safety.."];Quit[]];

outlabelsL={
	"\!\(\*SubsuperscriptBox[\(\[Theta]\), \(12\), \(l\)]\)",
	"\!\(\*SubsuperscriptBox[\(\[Theta]\), \(13\), \(l\)]\)",
	"\!\(\*SubsuperscriptBox[\(\[Theta]\), \(23\), \(l\)]\)",
	"\!\(\*SuperscriptBox[\(\[Delta]\), \(l\)]\)",
	"\!\(\*SubsuperscriptBox[\(\[CapitalDelta]m\), \(21\), \(2\)]\)",
	"\!\(\*SubsuperscriptBox[\(\[CapitalDelta]m\), \(31\), \(2\)]\)",
	"\!\(\*SubscriptBox[\(y\), \(e\)]\)",
	"\!\(\*SubscriptBox[\(y\), \(\[Mu]\)]\)",
	"\!\(\*SubscriptBox[\(y\), \(\[Tau]\)]\)"
};
outlabelsQ={
	"\!\(\*SubsuperscriptBox[\(\[Theta]\), \(12\), \(q\)]\)",
	"\!\(\*SubsuperscriptBox[\(\[Theta]\), \(13\), \(q\)]\)",
	"\!\(\*SubsuperscriptBox[\(\[Theta]\), \(23\), \(q\)]\)",
	"\!\(\*SuperscriptBox[\(\[Delta]\), \(q\)]\)",
	"\!\(\*SubscriptBox[\(y\), \(u\)]\)",
	"\!\(\*SubscriptBox[\(y\), \(c\)]\)",
	"\!\(\*SubscriptBox[\(y\), \(t\)]\)",
	"\!\(\*SubscriptBox[\(y\), \(d\)]\)",
	"\!\(\*SubscriptBox[\(y\), \(s\)]\)",
	"\!\(\*SubscriptBox[\(y\), \(b\)]\)"
};
outLabels=Join[outlabelsQ,outlabelsL];

findCredibilityLevel[grid_,level_]:=Module[{srtdata,cumsum,index},
	srtdata=grid//Flatten//Sort;
	cumsum=Accumulate[srtdata];
	(*overcover interval with Floor; Ceil would undercover*)
	index=Floor[BinarySearch[cumsum,(1-level)Last[cumsum]]];
	srtdata[[index]]
];


End[];
EndPackage[];
