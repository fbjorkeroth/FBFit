(* ::Package:: *)

BeginPackage["FBFit`Analysis`",{"FBFit`CalculateParameters`","FBFit`BestFitsAndErrors`"}];
Quiet[Needs["Combinatorica`"],General::compat];

FBImportFrom::usage="Sets the data directory to import from.";
FBExtractBestInput::usage="Extracts from the dataset the input parameters for the best fit point.";

FBPrintInput::usage="prints the best fit \!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\) and input parameters.";
FBPrintOutput::usage="prints the best fit physical parameters.";
FBPlotPulls::usage="plots a bar chart of the pulls of each output parameter.";

FBPhysicalParameterTable::usage="Outputs the physical parameters (Yukawa eigenvalues, mixing parameters) 
	for a given set of input data. Option exists to thin this data for quicker evaluation.";
FBCredibleInterval::usage="Calculates the N% credible interval for physical (Yukawa, mixing) parameters for a given input data set.";
FBPlotHistogram::usage="Plots histogram of physical output data.";

FBChopDataFraction::usage="Chops off the first Nth fraction of rows in a data set, where 0<N<1.";

Begin["Private`"];

(* ::Global variables:: *)

inLabels = Global`InLabels;

(* ::Function options:: *)

Options[FBPrintInput]={"Model"->"MSSM","TanB"->5.,"EtaB"->0.};
Options[FBPrintOutput]={"Model"->"MSSM","TanB"->5.,"EtaB"->0.,"Sector"->"All"};
Options[FBPlotPulls]={"Model"->"MSSM","TanB"->5.,"EtaB"->0.,"Sector"->"All"};

Options[FBPhysicalParameterTable]={"Thinning"->1};
Options[FBCredibleInterval]={"Thinning"->1,"SigmaSpan"->3,"PixelDensity"->100,"Plot"->True};
Options[FBPlotHistogram]={"Thinning"->1,"ImageSize"->300,"Bins"->30};

(* ::Public functions:: *)

FBImportFrom[runname_]:=Module[{dir},
	dir=FileNameJoin@{NotebookDirectory[],"data",runname};
	If[DirectoryQ[dir]==False,
		Print["FBImportFrom: data directory ",dir," not found! Please select another run. Quitting kernel for safety.."];
		Quit[]
	];
	SetDirectory[dir];
	Print["FBImportFrom: directory set to ",dir]
];

FBExtractBestInput[inputdata_]:=Module[{chisq,best},
	chisq=inputdata[[;;,-1]];
	best=Position[chisq,Min[chisq]][[1,1]];
	inputdata[[best,;;-2]]
];

FBPrintInput[t_,OptionsPattern[]]:=Module[{databestfit,dataerrors,titles,chisq,table},
	databestfit=FBGetDataBestFit[];
	dataerrors=FBGetDataErrors[];
	titles={"Parameter","Value"};
	chisq=FBChiSq[FBGetPulls[FBGetPhysicalParameters[t],databestfit,dataerrors]];	
(*	table=MatrixForm@Prepend[Transpose@{inLabels[[#]],t[[#]]},titles]&/@{isQuark,isLepton};*)
(*	Print[
		"\!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\): ",chisq,
		"\nInput: (quarks) ",table[[1]],"\t(leptons) ",table[[2]]
	]*)
	table=MatrixForm@Prepend[Transpose@{inLabels,t},titles];
	Print[
		"\!\(\*SuperscriptBox[\(\[Chi]\), \(2\)]\): ",chisq,
		"\nInput: ",table
	]
];

FBPrintOutput[t_,OptionsPattern[]]:=Module[{databestfit,calc,titles,table,x},
	databestfit=FBGetDataBestFit[];
	calc=FBGetPhysicalParameters[t];
	titles={"Parameter","Data","Model"};
	table=MatrixForm@Prepend[Transpose@{outLabels[OptionValue["Sector"]],databestfit,calc},titles];
	x=ScientificForm[table,4,ExponentFunction->(If[-2<#<3,Null,#]&)];
	Print["Output: ",x]
];

FBPlotPulls[t_,OptionsPattern[]]:=Module[{databestfit,dataerrors,calc,pulls},
	databestfit=FBGetDataBestFit[];
	dataerrors=FBGetDataErrors[];
	calc=FBGetPhysicalParameters[t];
	pulls=FBGetPulls[calc,databestfit,dataerrors];
	Print@BarChart[pulls,ChartLabels->outLabels[OptionValue["Sector"]],AxesLabel->"Pull",ImageSize->Large,AspectRatio->1/2,BaseStyle->FontSize->12]
];

FBPhysicalParameterTable[inputdata_,variable_,OptionsPattern[]]:=Module[{dataThinned,n=variable},
	dataThinned=Take[inputdata,{1,Length[inputdata],OptionValue["Thinning"]}];
	FBGetPhysicalParameters[#][[n]]&/@dataThinned[[;;,;;-2]]
];


FBCredibleInterval[inputdata_,n_,CIlevel_,OptionsPattern[]]:=Module[{l,mu,std,pdf,spacing,grid,t,plot,ci},
	l=FBPhysicalParameterTable[inputdata,n,"Thinning"->OptionValue["Thinning"]];
	{mu,std,pdf}=#[l]&/@{Mean,StandardDeviation,SmoothKernelDistribution};
	spacing=10std/OptionValue["PixelDensity"];
	grid=Table[PDF[pdf,{mu-OptionValue["SigmaSpan"] std+spacing*i}],{i,OptionValue["PixelDensity"]}];
	t=findCredibilityLevel[grid,CIlevel];

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

FBPlotHistogram[inputdata_,n_,OptionsPattern[]]:=Module[{bf,err,p,bins,h,maxH,line},
	bf=FBGetDataBestFit[][[n]];
	err=FBGetDataErrors[][[n]];
	p=FBPhysicalParameterTable[inputdata,n,"Thinning"->OptionValue["Thinning"]];
	{bins, h} = HistogramList[p,OptionValue["Bins"]];
	maxH=Max[h];
	line[x_]:=Line[{{x,0},{x,maxH+2}}];
	Histogram[
		p,{bins},h&,
		Epilog->{line[bf],{Dashed,line[bf+err]},{Dashed,line[bf-err]}}
	]
];

FBChopDataFraction[inputdata_,bottomfraction_]:=Module[{n},
	n=Round[(1-bottomfraction)Length[inputdata]];
	Take[inputdata,-n]
];


(* ::Internal functions:: *)

printBadModel[]:=Module[{},Print["Model not supported. Quitting kernel for safety.."];Quit[]];

outLabels[sec_]:=Module[{q,l},
	q={
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
	l={
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
	Switch[sec,
		"Q",q,
		"L",l,
		sec,Join[q,l]
	]
];

findCredibilityLevel[grid_,level_]:=Module[{srtdata,cumsum,index},
	srtdata=grid//Flatten//Sort;
	cumsum=Accumulate[srtdata];
	(*overcover interval with Floor; Ceil would undercover*)
	index=Floor[BinarySearch[cumsum,(1-level)Last[cumsum]]];
	srtdata[[index]]
];


End[];
EndPackage[];
