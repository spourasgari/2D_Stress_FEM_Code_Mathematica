(* ::Package:: *)

ClearAll["Global`*"]

tool=2;
arzfull=1;
arz=0.5;
hnum=4;
vnumtot=8;
vnum=vnumtot/2;
radius=0.25;
EM=200*10^9;
v=0.3;
h=0.05*10^(-3);
P=1000;
affectedwidth=0.5;
Ftot=P*h*affectedwidth;
(*hnum=Input["Enter the number of nodes in x direction on the upper edge"]
vnum=Input["Enter the number of nodes in x direction on the upper edge"]*)

hint=tool/(hnum-1);
vint=arz/(vnum-1);

EleCountV=(vnum-1)*2;
EleCountH=(hnum-1)*2;

EleCount=(hnum-1)*(vnum-1)*2;

NodalConnect=ConstantArray[0,{EleCount,4}];
grid=ConstantArray[0,{EleCount,3}];
For[EleNum=1, EleNum<=EleCount,  EleNum++,
		rownum=If[Mod[EleNum,EleCountV]==0,EleCountV/2,Ceiling[Mod[EleNum,EleCountV]/2]];
		colnum=Ceiling[EleNum/EleCountV];
		grid[[EleNum,1]]=EleNum;
		grid[[EleNum,2]]=rownum;
		grid[[EleNum,3]]=colnum;
		
		NodalConnect[[EleNum,1]]=EleNum;
		NodalConnect[[EleNum,2]]=rownum+(colnum-1)*vnum+1;
		NodalConnect[[EleNum,3]]=rownum+(colnum)*vnum;
		NodalConnect[[EleNum,4]]=If[Mod[EleNum,2]==1,rownum+(colnum-1)*vnum, rownum+(colnum-1)*vnum+1+vnum]
]

NodeCount=(rownum+(colnum-1)*vnum+1+vnum)/.{rownum->grid[[EleCount,2]],colnum->grid[[EleCount,3]]};

NodeCoord=ConstantArray[0,{NodeCount,3}];

NodeNum=1;
x=0;
	For[i=1,i<=hnum,i++,
		For[j=1,j<=vnum,j++,
		y=arz-(j-1)*vint;
		x=(i-1)*hint;
			
			NodeCoord[[NodeNum,1]]=NodeNum;
			NodeCoord[[NodeNum,2]]=x;
			NodeCoord[[NodeNum,3]]=y;
			
			NodeNum++
		]
	];


(*Modifying the mesh around hole*)
(*radialdistance=Sqrt((NodeCoord[[NodeNum,2]]-tool/2)^2+(NodeCoord[[NodeNum,3]])^2)
For[NodeNum=1,NodeNum\[LessEqual]NodeCount,NodeNum++,
	If[radialdistance\[LessEqual]radius,*)


SF=ConstantArray[0,{EleCount,8}];

(* Matrix SF: Element, Node1, Node2, Node3, SF1, SF2, SF3, A_of_element *)

For[i=1,i<=EleCount,i++,
	For[j=1,j<=4,j++,
		SF[[i,j]]=NodalConnect[[i,j]]
		];
	SF[[i,8]]=-0.5*Det[
					{{1,NodeCoord[[NodalConnect[[i,2]],2]],NodeCoord[[NodalConnect[[i,2]],3]]},
					{1,NodeCoord[[NodalConnect[[i,3]],2]],NodeCoord[[NodalConnect[[i,3]],3]]},
					{1,NodeCoord[[NodalConnect[[i,4]],2]],NodeCoord[[NodalConnect[[i,4]],3]]}}
					];

	SF[[i,5]]=-0.5*Det[
					{{1,r,s},
					{1,NodeCoord[[NodalConnect[[i,3]],2]],NodeCoord[[NodalConnect[[i,3]],3]]},
					{1,NodeCoord[[NodalConnect[[i,4]],2]],NodeCoord[[NodalConnect[[i,4]],3]]}}
					]/SF[[i,8]];

	SF[[i,6]]=-0.5*Det[
					{{1,NodeCoord[[NodalConnect[[i,2]],2]],NodeCoord[[NodalConnect[[i,2]],3]]},
					{1,r,s},
					{1,NodeCoord[[NodalConnect[[i,4]],2]],NodeCoord[[NodalConnect[[i,4]],3]]}}
					]/SF[[i,8]];

	SF[[i,7]]=-0.5*Det[
					{{1,NodeCoord[[NodalConnect[[i,2]],2]],NodeCoord[[NodalConnect[[i,2]],3]]},
					{1,NodeCoord[[NodalConnect[[i,3]],2]],NodeCoord[[NodalConnect[[i,3]],3]]},
					{1,r,s}}
					]/SF[[i,8]]
	]


Kg=ConstantArray[0,{2*NodeCount,2*NodeCount}];
EMat={{EM/(1-v^2),v*EM/(1-v^2),0},{v*EM/(1-v^2),EM/(1-v^2),0},{0,0,0.5*EM/(1+v)}};


For[EleNum=1,EleNum<=EleCount,EleNum++,
	B={{D[SF[[EleNum,5]],r],0,D[SF[[EleNum,6]],r],0,D[SF[[EleNum,7]],r],0},
		{0,D[SF[[EleNum,5]],s],0,D[SF[[EleNum,6]],s],0,D[SF[[EleNum,7]],s]},
		{D[SF[[EleNum,5]],s],D[SF[[EleNum,5]],r],D[SF[[EleNum,6]],s],D[SF[[EleNum,6]],r],D[SF[[EleNum,7]],s],D[SF[[EleNum,7]],r]}};

		Ke=Transpose[B].EMat.B*SF[[EleNum,8]];

		Kg[[NodalConnect[[EleNum,2]]*2-1,NodalConnect[[EleNum,2]]*2-1]]+=Ke[[1,1]];
		Kg[[NodalConnect[[EleNum,2]]*2-1,NodalConnect[[EleNum,2]]*2]]+=Ke[[1,2]];
		Kg[[NodalConnect[[EleNum,2]]*2,NodalConnect[[EleNum,2]]*2-1]]+=Ke[[2,1]];
		Kg[[NodalConnect[[EleNum,2]]*2,NodalConnect[[EleNum,2]]*2]]+=Ke[[2,2]];

		Kg[[NodalConnect[[EleNum,3]]*2-1,NodalConnect[[EleNum,3]]*2-1]]+=Ke[[3,3]];
		Kg[[NodalConnect[[EleNum,3]]*2-1,NodalConnect[[EleNum,3]]*2]]+=Ke[[3,4]];
		Kg[[NodalConnect[[EleNum,3]]*2,NodalConnect[[EleNum,3]]*2-1]]+=Ke[[4,3]];
		Kg[[NodalConnect[[EleNum,3]]*2,NodalConnect[[EleNum,3]]*2]]+=Ke[[4,4]];

		Kg[[NodalConnect[[EleNum,4]]*2-1,NodalConnect[[EleNum,4]]*2-1]]+=Ke[[5,5]];
		Kg[[NodalConnect[[EleNum,4]]*2-1,NodalConnect[[EleNum,4]]*2]]+=Ke[[5,6]];
		Kg[[NodalConnect[[EleNum,4]]*2,NodalConnect[[EleNum,4]]*2-1]]+=Ke[[6,5]];
		Kg[[NodalConnect[[EleNum,4]]*2,NodalConnect[[EleNum,4]]*2]]+=Ke[[6,6]];


		Kg[[NodalConnect[[EleNum,2]]*2-1,NodalConnect[[EleNum,3]]*2-1]]+=Ke[[1,3]];
		Kg[[NodalConnect[[EleNum,2]]*2-1,NodalConnect[[EleNum,3]]*2]]+=Ke[[1,4]];
		Kg[[NodalConnect[[EleNum,2]]*2,NodalConnect[[EleNum,3]]*2-1]]+=Ke[[2,3]];
		Kg[[NodalConnect[[EleNum,2]]*2,NodalConnect[[EleNum,3]]*2]]+=Ke[[2,4]];

		Kg[[NodalConnect[[EleNum,2]]*2-1,NodalConnect[[EleNum,4]]*2-1]]+=Ke[[1,5]];
		Kg[[NodalConnect[[EleNum,2]]*2-1,NodalConnect[[EleNum,4]]*2]]+=Ke[[1,6]];
		Kg[[NodalConnect[[EleNum,2]]*2,NodalConnect[[EleNum,4]]*2-1]]+=Ke[[2,5]];
		Kg[[NodalConnect[[EleNum,2]]*2,NodalConnect[[EleNum,4]]*2]]+=Ke[[2,6]];

		Kg[[NodalConnect[[EleNum,3]]*2-1,NodalConnect[[EleNum,2]]*2-1]]+=Ke[[3,1]];
		Kg[[NodalConnect[[EleNum,3]]*2-1,NodalConnect[[EleNum,2]]*2]]+=Ke[[3,2]];
		Kg[[NodalConnect[[EleNum,3]]*2,NodalConnect[[EleNum,2]]*2-1]]+=Ke[[4,1]];
		Kg[[NodalConnect[[EleNum,3]]*2,NodalConnect[[EleNum,2]]*2]]+=Ke[[4,2]];

		Kg[[NodalConnect[[EleNum,3]]*2-1,NodalConnect[[EleNum,4]]*2-1]]+=Ke[[3,5]];
		Kg[[NodalConnect[[EleNum,3]]*2-1,NodalConnect[[EleNum,4]]*2]]+=Ke[[3,6]];
		Kg[[NodalConnect[[EleNum,3]]*2,NodalConnect[[EleNum,4]]*2-1]]+=Ke[[4,5]];
		Kg[[NodalConnect[[EleNum,3]]*2,NodalConnect[[EleNum,4]]*2]]+=Ke[[4,6]];

		Kg[[NodalConnect[[EleNum,4]]*2-1,NodalConnect[[EleNum,2]]*2-1]]+=Ke[[5,1]];
		Kg[[NodalConnect[[EleNum,4]]*2-1,NodalConnect[[EleNum,2]]*2]]+=Ke[[5,2]];
		Kg[[NodalConnect[[EleNum,4]]*2,NodalConnect[[EleNum,2]]*2-1]]+=Ke[[6,1]];
		Kg[[NodalConnect[[EleNum,4]]*2,NodalConnect[[EleNum,2]]*2]]+=Ke[[6,2]];

		Kg[[NodalConnect[[EleNum,4]]*2-1,NodalConnect[[EleNum,3]]*2-1]]+=Ke[[5,3]];
		Kg[[NodalConnect[[EleNum,4]]*2-1,NodalConnect[[EleNum,3]]*2]]+=Ke[[5,4]];
		Kg[[NodalConnect[[EleNum,4]]*2,NodalConnect[[EleNum,3]]*2-1]]+=Ke[[6,3]];
		Kg[[NodalConnect[[EleNum,4]]*2,NodalConnect[[EleNum,3]]*2]]+=Ke[[6,4]]

]
KgFinal=0.5*h*Kg;





ClearAll[fsupx,fsupy,fsupport,fsym,fsymmetry]
(*Distribiution of total force on affected nodes (rounded up!) *)
affnodecount=Ceiling[(affectedwidth/arz)*vnum];
FNode=Ftot/affnodecount;

(*fsupport=ConstantArray[0,{vnum,2}];*)

(*fsym=ConstantArray[0,{hnum,1}]; *)  (*Number of node, unknown vertical force*)

nodesonsym=ConstantArray[0,{hnum,1}];

j=1;
For[i=1,i<=NodeCount,i++,
	If[NodeCoord[[i,3]]==0,nodesonsym[[j]]=i;j++]
]


FMat=ConstantArray[0,{2*NodeCount,1}];

(*Force on nodes*)

(*near force*)
For[i=1,i<=NodeCount,i++,
		If[NodeCoord[[i,2]]==tool && NodeCoord[[i,3]]<=(affectedwidth/2+vint/2),FMat[[2*i-1]]={FNode}];
]

(* on symmetry line*)
j=1;
Do[
	FMat[[2*i]]={fsym[i]};
	j++, {i,nodesonsym}
]

(*near support*)
For[i=1,i<=NodeCount,i++,
		If[NodeCoord[[i,2]]==0,FMat[[2*i-1]]={fsupx[i]};FMat[[2*i]]={fsupy[i]}]
]





UMat=Inverse[KgFinal].FMat;
constraints=ConstantArray[0,{2*vnum+(hnum-1)}];
j=1;
Do[constraints[[j]]=UMat[[i*2]]==0; j++, {i,nodesonsym}];
Do[constraints[[j]]=UMat[[i*2]]==0; j++, {i,vnum-1}];
Do[constraints[[j]]=UMat[[i*2-1]]==0; j++, {i,vnum}];


unknownf=ConstantArray[0,{2*vnum+(hnum-1)}];
j=1;

Do[unknownf[[j]]=fsupx[i]; j++, {i,vnum}];
Do[unknownf[[j]]=fsupy[i]; j++, {i,vnum}];
j-=1;
Do[unknownf[[j]]=fsym[i]; j++, {i,nodesonsym}];





unknownff=unknownf/.{fsym[vnum]->fsupy[4]}


ForceResult=Solve[constraints,unknownff];
DispResult=((UMat/.{ForceResult})[[1,1]]);





EpsMat=ConstantArray[0,{EleCount,3}];
StressMat=ConstantArray[0,{EleCount,3}];
vonMises=ConstantArray[0,{EleCount,1}];

For[EleNum=1,EleNum<=EleCount,EleNum++,
	B={{D[SF[[EleNum,5]],r],0,D[SF[[EleNum,6]],r],0,D[SF[[EleNum,7]],r],0},
		{0,D[SF[[EleNum,5]],s],0,D[SF[[EleNum,6]],s],0,D[SF[[EleNum,7]],s]},
		{D[SF[[EleNum,5]],s],D[SF[[EleNum,5]],r],D[SF[[EleNum,6]],s],D[SF[[EleNum,6]],r],D[SF[[EleNum,7]],s],D[SF[[EleNum,7]],r]}};
		
			Ue={
			DispResult[[NodalConnect[[EleNum,2]]*2-1]],
			DispResult[[NodalConnect[[EleNum,2]]*2]],
			DispResult[[NodalConnect[[EleNum,3]]*2-1]],
			DispResult[[NodalConnect[[EleNum,3]]*2]],
			DispResult[[NodalConnect[[EleNum,4]]*2-1]],
			DispResult[[NodalConnect[[EleNum,4]]*2]]};

			epsE=B.Ue;
(*Global Epsilon Matrix= EleCount*3: exx, eyy, 2exy*)
			EpsMat[[EleNum,1]]=epsE[[1,1]];
			EpsMat[[EleNum,2]]=epsE[[2,1]];
			EpsMat[[EleNum,3]]=epsE[[3,1]];

			strssE=EMat.B.Ue;
			StressMat[[EleNum,1]]=strssE[[1,1]];
			StressMat[[EleNum,2]]=strssE[[2,1]];
			StressMat[[EleNum,3]]=strssE[[3,1]];

			vonMises[[EleNum,1]]=Sqrt[StressMat[[EleNum,1]]^2+StressMat[[EleNum,2]]^2-StressMat[[EleNum,1]]*StressMat[[EleNum,2]]+3*StressMat[[EleNum,3]]^2]
]





StressCoord=ConstantArray[0,{EleCount,3}];
For[EleNum=1,EleNum<=EleCount,EleNum++,
	StressCoord[[EleNum,3]]=vonMises[[EleNum,1]];
	StressCoord[[EleNum,1]]=(NodeCoord[[NodalConnect[[EleNum,2]],2]]+NodeCoord[[NodalConnect[[EleNum,3]],2]]+NodeCoord[[NodalConnect[[EleNum,4]],2]])/3;
	StressCoord[[EleNum,2]]=(NodeCoord[[NodalConnect[[EleNum,2]],3]]+NodeCoord[[NodalConnect[[EleNum,3]],3]]+NodeCoord[[NodalConnect[[EleNum,4]],3]])/3
]


ListDensityPlot[StressCoord]
