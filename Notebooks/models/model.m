Yu={{Cu11, Cu12 Exp[I z12], 0},{0, Cu22, 0},{0, 0, Cu33}}//ConjugateTranspose;
Yd={{Cd11, Cd12, 0},{0, Cd22, Cd23},{0, 0, Cd33}}//ConjugateTranspose;
Ye=DiagonalMatrix@({4.90856087*^-4,0.103622931,1.76167}/174);
Mnu=(
	Mu1{{0, 0, 0},{0, 1, 1},{0, 1, 1}}
	+ Mu2 E^(I zMu2){{1, 4, 2},{4, 16, 8},{2, 8, 4}}
	+ Mu3 E^(I zMu3) {{0, 0, 0},{0, 0, 0},{0, 0, 1}}
)/.{Mu1->0.035,Mu2->0.002,Mu3->0.0015,zMu2->4\[Pi]/5,zMu3->0}//ConjugateTranspose;


InputVariables={Cu11,Cu12,Cu22,Cu33,Cd11,Cd12,Cd22,Cd23,Cd33,z12};
InLabels=ToString/@InputVariables;

IsReal=Range[1,9];
IsPhase={10};
IsQuark=Range[Length[inputVariables]];
IsLepton={};

StartBounds={
	{0.00001,0.0001},{0.0001,0.001},{0.001,0.01},{0.5,1},
	{0.00001,0.0001},{0.0001,0.001},{0.0001,0.001},{0.001,0.01},{0.01,0.1},
	{0,2\[Pi]}
};