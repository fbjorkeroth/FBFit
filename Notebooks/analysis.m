(* ::Package:: *)

(* ::Title:: *)
(*MCMC analysis*)


(* ::Section:: *)
(*Load packages and data*)


Get["FBFit`"];

FBLoadModel["models/C8.m"];
FBImportFrom["test"];

likelihoodTable=Import["rundata.txt","Table"];
acceptanceList=Import["acceptance.log","List"];
sigmaList=Import["sigma.log","List"];
options=Import["variables.mx"];

FBSetOptions@@options;
FBLoadBestFitsAndErrors[];


(* ::Section:: *)
(*Analysis*)


(* ::Subsection:: *)
(*Adjust data table*)


likelihoodTable=FBChopDataFraction[likelihoodTable,0.25];


(* ::Subsection:: *)
(*Evolution of \[Chi]^2*)


chisqList=likelihoodTable[[;;,-1]];
ListPlot[chisqList]


(* ::Subsection:: *)
(*Best fit point*)


thetaBest=FBExtractBestInput[likelihoodTable];

FBPrintInput[thetaBest];
FBPrintOutput[thetaBest];
FBPlotPulls[thetaBest];


(* ::Subsubsection:: *)
(*Reconstruct matrices*)


yun=Yu/.Thread[InputVariables->thetaBest];
yun//MatrixForm
ydn=Yd/.Thread[InputVariables->thetaBest];
ydn//MatrixForm


(*Quiet@Needs["MixingParameterTools`MPT3x3`"];*)
CKMParameters[yun,DiagonalMatrix@{1,2,3}][[1]]
CKMParameters[ydn,DiagonalMatrix@{1,2,3}][[1]]
CKMParameters[yun,ydn][[1]]


(* ::Subsection:: *)
(*Evolution of links (\[Alpha], \[Sigma])*)


range=1000;
ma=ListLogLinearPlot[Table[Mean[acceptanceList[[;;n]]],{n,range}],GridLines->{{},{0.3,0.5}},PlotRange->All,Joined->True];
si=ListLogLinearPlot[40sigmaList,PlotRange->{{1,range},All},PlotStyle->Orange];
ac=ListLogLinearPlot[acceptanceList,PlotRange->{{1,range},All}];
Show[{ma,si}]


(* ::Subsection:: *)
(*Output parameter distributions*)


(* ::Subsubsection:: *)
(*Histogram*)


FBPlotHistogram[likelihoodTable,7,"Thinning"->20,"Bins"->50]


(* ::Subsubsection:: *)
(*Credible interval*)


(*FBCredibleInterval[likelihoodTable,#,0.95,"Thinning"\[Rule]100]&/@Range[1,10];*)
