/* kmetis.c */
#define	Mc_Global_Partition Mc_Global_Partition__

/* mmetis.c */

/* gkmetis.c */

/* match.c */
#define	Mc_GlobalMatch_Balance Mc_GlobalMatch_Balance__

/* coarsen.c */
#define	Mc_Global_CreateCoarseGraph Mc_Global_CreateCoarseGraph__

/* initpart.c */
#define	Mc_InitPartition_RB Mc_InitPartition_RB__
#define	Mc_KeepPart Mc_KeepPart__

/* kwayrefine.c */
#define	Mc_ProjectPartition Mc_ProjectPartition__
#define	Mc_ComputePartitionParams Mc_ComputePartitionParams__

/* kwayfm.c */
#define	Mc_KWayFM Mc_KWayFM__

/* kwaybalance.c */
#define	Mc_KWayBalance Mc_KWayBalance__

/* remap.c */
#define	ParallelReMapGraph ParallelReMapGraph__
#define	ParallelTotalVReMap ParallelTotalVReMap__
#define	SimilarTpwgts SimilarTpwgts__

/* move.c */
#define	Mc_MoveGraph Mc_MoveGraph__
#define	CheckMGraph CheckMGraph__
#define	ProjectInfoBack ProjectInfoBack__
#define	FindVtxPerm FindVtxPerm__

/* memory.c */
#define	PreAllocateMemory PreAllocateMemory__
#define	FreeWSpace FreeWSpace__
#define	FreeCtrl FreeCtrl__
#define	CreateGraph CreateGraph__
#define	InitGraph InitGraph__
#define	FreeGraph FreeGraph__
#define	FreeInitialGraphAndRemap FreeInitialGraphAndRemap__


/************************/
/* Adaptive subroutines */
/************************/
/* ametis.c */
#define	Adaptive_Partition Adaptive_Partition__

/* rmetis.c */

/* lmatch.c */
#define	Mc_LocalMatch_HEM Mc_LocalMatch_HEM__
#define	Mc_Local_CreateCoarseGraph Mc_Local_CreateCoarseGraph__

/* wave.c */
#define	WavefrontDiffusion WavefrontDiffusion__

/* balancemylink.c */
#define	BalanceMyLink BalanceMyLink__

/* redomylink.c */
#define	RedoMyLink RedoMyLink__

/* initbalance.c */
#define	Balance_Partition Balance_Partition__
#define	Mc_AssembleAdaptiveGraph Mc_AssembleAdaptiveGraph__

/* mdiffusion.c */
#define	Mc_Diffusion Mc_Diffusion__
#define	ExtractGraph ExtractGraph__

/* diffutil.c */
#define	SetUpConnectGraph SetUpConnectGraph__
#define	Mc_ComputeMoveStatistics Mc_ComputeMoveStatistics__
#define	Mc_ComputeSerialTotalV Mc_ComputeSerialTotalV__
#define	ComputeLoad ComputeLoad__
#define	ConjGrad2 ConjGrad2__
#define	mvMult2 mvMult2__
#define	ComputeTransferVector ComputeTransferVector__
#define	ComputeSerialEdgeCut ComputeSerialEdgeCut__
#define	ComputeSerialTotalV ComputeSerialTotalV__

/* akwayfm.c */
#define	Mc_KWayAdaptiveRefine Mc_KWayAdaptiveRefine__

/* selectq.c */
#define	Mc_DynamicSelectQueue Mc_DynamicSelectQueue__
#define	Mc_HashVwgts Mc_HashVwgts__
#define	Mc_HashVRank Mc_HashVRank__

/* csrmatch.c */
#define	CSR_Match_SHEM CSR_Match_SHEM__

/* serial.c */
#define	Mc_SerialKWayAdaptRefine Mc_SerialKWayAdaptRefine__
#define	Mc_ComputeSerialPartitionParams Mc_ComputeSerialPartitionParams__
#define	AreAllHVwgtsBelow AreAllHVwgtsBelow__
#define	ComputeHKWayLoadImbalance ComputeHKWayLoadImbalance__
#define	SerialRemap SerialRemap__
#define	SSMIncKeyCmp SSMIncKeyCmp__
#define	Mc_Serial_FM_2WayRefine Mc_Serial_FM_2WayRefine__
#define	Serial_SelectQueue Serial_SelectQueue__
#define	Serial_BetterBalance Serial_BetterBalance__
#define	Serial_Compute2WayHLoadImbalance Serial_Compute2WayHLoadImbalance__
#define	Mc_Serial_Balance2Way Mc_Serial_Balance2Way__
#define	Mc_Serial_Init2WayBalance Mc_Serial_Init2WayBalance__
#define	Serial_SelectQueueOneWay Serial_SelectQueueOneWay__
#define	Mc_Serial_Compute2WayPartitionParams Mc_Serial_Compute2WayPartitionParams__
#define	Serial_AreAnyVwgtsBelow Serial_AreAnyVwgtsBelow__

/* weird.c */
#define	PartitionSmallGraph PartitionSmallGraph__
#define	CheckInputs CheckInputs__


/****************************/
/* Mesh to Dual subroutines */
/****************************/
/* mesh.c */
/* msetup.c */
#define	SetUpMesh SetUpMesh__
#define	CreateMesh CreateMesh__
#define	InitMesh InitMesh__


/************************/
/* Ordering subroutines */
/************************/
/* ometis.c */
/* pspases.c */
#define	AssembleEntireGraph AssembleEntireGraph__

/* node_refine.c */
#define	ComputeNodePartitionParams0 ComputeNodePartitionParams0__
#define	ComputeNodePartitionParams ComputeNodePartitionParams__
#define	KWayNodeRefine0 KWayNodeRefine0__
#define	KWayNodeRefine KWayNodeRefine__
#define	KWayNodeRefine2 KWayNodeRefine2__
#define	PrintNodeBalanceInfo PrintNodeBalanceInfo__

/* initmsection.c */
#define	InitMultisection InitMultisection__
#define	AssembleMultisectedGraph AssembleMultisectedGraph__

/* order.c */
#define	MultilevelOrder MultilevelOrder__
#define	LabelSeparators LabelSeparators__
#define	CompactGraph CompactGraph__
#define	LocalOrder LocalOrder__
#define	LocalNDOrder LocalNDOrder__
#define	Order_Partition Order_Partition__

/* xyzpart.c */
#define	Coordinate_Partition Coordinate_Partition__
#define	PartSort PartSort__

/***********************/
/* Utility subroutines */
/***********************/
/* fpqueue.c */
#define	FPQueueInit FPQueueInit__
#define	FPQueueReset FPQueueReset__
#define	FPQueueFree FPQueueFree__
#define	FPQueueGetSize FPQueueGetSize__
#define	FPQueueInsert FPQueueInsert__
#define	FPQueueDelete FPQueueDelete__
#define	FPQueueUpdate FPQueueUpdate__
#define	FPQueueUpdateUp FPQueueUpdateUp__
#define	FPQueueGetMax FPQueueGetMax__
#define	FPQueueSeeMaxVtx FPQueueSeeMaxVtx__
#define	FPQueueSeeMaxGain FPQueueSeeMaxGain__
#define	FPQueueGetKey FPQueueGetKey__
#define	FPQueueGetQSize FPQueueGetQSize__
#define	CheckHeapFloat CheckHeapFloat__

/* stat.c */
#define	Mc_ComputeSerialBalance Mc_ComputeSerialBalance__
#define	Mc_ComputeParallelBalance Mc_ComputeParallelBalance__
#define	Mc_PrintThrottleMatrix Mc_PrintThrottleMatrix__
#define	Mc_ComputeRefineStats Mc_ComputeRefineStats__

/* debug.c */
#define	PrintVector PrintVector__
#define	PrintVector2 PrintVector2__
#define	PrintPairs PrintPairs__
#define	PrintGraph PrintGraph__
#define	PrintGraph2 PrintGraph2__
#define	PrintSetUpInfo PrintSetUpInfo__
#define	PrintTransferedGraphs PrintTransferedGraphs__
#define	WriteMetisGraph WriteMetisGraph__

/* comm.c */
#define	CommInterfaceData CommInterfaceData__
#define	CommChangedInterfaceData CommChangedInterfaceData__
#define	GlobalSEMax GlobalSEMax__
#define	GlobalSEMaxDouble GlobalSEMaxDouble__
#define	GlobalSEMin GlobalSEMin__
#define	GlobalSESum GlobalSESum__
#define	GlobalSEMaxFloat GlobalSEMaxFloat__
#define	GlobalSEMinFloat GlobalSEMinFloat__
#define	GlobalSESumFloat GlobalSESumFloat__

/* util.c */
#define	errexit errexit__
#define	myprintf myprintf__
#define	rprintf rprintf__
#define	imalloc imalloc__
#define	idxmalloc idxmalloc__
#define	fmalloc fmalloc__
#define	ismalloc ismalloc__
#define	idxsmalloc idxsmalloc__
#define	GKmalloc GKmalloc__
#define	GKfree GKfree__
#define	iset iset__
#define	idxset idxset__
#define	idxamax idxamax__
#define	idxamin idxamin__
#define	idxasum idxasum__
#define	snorm2 snorm2__
#define	sdot sdot__
#define	saxpy saxpy__
#define	ikeyvalsort_org ikeyvalsort_org__
#define	IncKeyValueCmp IncKeyValueCmp__
#define	dkeyvalsort dkeyvalsort__
#define	DecKeyValueCmp DecKeyValueCmp__
#define	BSearch BSearch__
#define	RandomPermute RandomPermute__
#define	FastRandomPermute FastRandomPermute__
#define	ispow2 ispow2__
#define	log2Int log2Int__
#define	BucketSortKeysDec BucketSortKeysDec__
#define	sset sset__
#define	iamax iamax__
#define	idxamax_strd idxamax_strd__
#define	idxamin_strd idxamin_strd__
#define	samax_strd samax_strd__
#define	sfamax sfamax__
#define	samin_strd samin_strd__
#define	idxavg idxavg__
#define	savg savg__
#define	samax samax__
#define	sfavg sfavg__
#define	samax2 samax2__
#define	samin samin__
#define	idxsum idxsum__
#define	idxsum_strd idxsum_strd__
#define	idxadd idxadd__
#define	ssum ssum__
#define	ssum_strd ssum_strd__
#define	sscale sscale__
#define	saneg saneg__
#define	BetterVBalance BetterVBalance__
#define	IsHBalanceBetterTT IsHBalanceBetterTT__
#define	IsHBalanceBetterFT IsHBalanceBetterFT__
#define	myvalkeycompare myvalkeycompare__
#define	imyvalkeycompare imyvalkeycompare__
#define	fsmalloc fsmalloc__
#define	saxpy2 saxpy2__
#define	GetThreeMax GetThreeMax__

/* qsort_special.c */
#define	iidxsort iidxsort__
#define	iintsort iintsort__
#define	ikeysort ikeysort__
#define	ikeyvalsort ikeyvalsort__

/* grsetup.c */
#define	Mc_SetUpGraph Mc_SetUpGraph__
#define	SetUpCtrl SetUpCtrl__
#define	ChangeNumbering ChangeNumbering__
#define	ChangeNumberingMesh ChangeNumberingMesh__
#define	GraphRandomPermute GraphRandomPermute__
#define	ComputeMoveStatistics ComputeMoveStatistics__

/* timer.c */
#define	InitTimers InitTimers__
#define	PrintTimingInfo PrintTimingInfo__
#define	PrintTimer PrintTimer__

/* setup.c */
#define	SetUp SetUp__
#define	Home_PE Home_PE__


