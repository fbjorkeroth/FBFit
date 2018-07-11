(* ::Title:: *)
(*MCMC analysis*)

(* ::Section:: *)
(*Load packages and data*)

Get["FBFit`"];

FBLoadModel["models/model.m"];

FBImportFrom["test"];

likelihoodTable=Import["rundata.txt","Table"];
acceptanceList=Import["acceptance.log","List"];
options=Import["variables.mx"];

FBSetOptions@@options;
FBLoadBestFitsAndErrors[];

(* ::Section:: *)
(*Plot best fit point*)

thetaBest=FBExtractBestInput[likelihoodTable];

FBPrintInput[thetaBest];
FBPrintOutput[thetaBest];
FBPlotPulls[thetaBest];

(* ::Section:: *)
(*Analytics*)

(* ::Subsection:: *)
(*Evolution of \[Alpha] (acceptance)*)

chisqList=likelihoodTable[[;;,-1]];

Histogram[chisqList]
chisqHist=Block[{ch,div=100,index,me},
	ch=chisqList/Mean[chisqList];
	index=Table[i,{i,Length[ch]}];
	Transpose[Mean/@Partition[#,div]&/@{index,ch}]
];
ListPlot[chisqHist,Joined->True]

ListPlot[MovingAverage[acceptanceList,100]]


(* ::Subsection:: *)
(*Evolution of mean \[Sigma]*)


checkt=List[];
Do[
	\[Sigma]mean=Table[Mean[acceptanceList[[;;n]]],{n,NNN}]//AbsoluteTiming;
	AppendTo[checkt,{NNN,\[Sigma]mean[[1]]}],
	{NNN,100,10000,100}
];
ListPlot[checkt]
(*ListLogLinearPlot[\[Sigma]mean[[2]],Joined->True,PlotRange->All]*)


(* ::Section:: *)
(*Reconstruct Yukawa matrices*)


Yu/.Thread[inputVariables->thetaBest]//Eigenvalues//Abs;
Yu/.Thread[inputVariables->thetaBest]//Abs//MatrixForm

Yd/.Thread[inputVariables->thetaBest]//Eigenvalues//Abs;
Yd/.Thread[inputVariables->thetaBest]//Abs//MatrixForm


(* ::Section:: *)
(*Bayesian analysis*)


FBPlotHistogram[likelihoodTable,1,"Thinning"->1,"Bins"->20]
(*FBCredibleInterval[likelihoodTable,#,0.95,"Thinning"\[Rule]100]&/@Range[1,10];*)
