/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ometis.c
 *
 * This is the entry point of parallel ordering
 *
 * Started 10/19/96
 * George
 *
 * $Id: ometis.c,v 1.4 2003/07/25 04:01:04 karypis Exp $
 *
 */

#include <parmetislib.h>




/***********************************************************************************
* This function is the entry point of the parallel ordering algorithm.
* This function assumes that the graph is already nice partitioned among the 
* processors and then proceeds to perform recursive bisection.
************************************************************************************/
void ParMETIS_V3_NodeND(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag,
              int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm)
{
  int i, j;
  int ltvwgts[MAXNCON];
  int nparts, npes, mype, wgtflag = 0, seed;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph, *mgraph;
  idxtype *morder;
  int minnvtxs;
  int dbglvl_original;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  nparts = 1*npes;

  if (!ispow2(npes)) {
    if (mype == 0)
      printf("Error: The number of processors must be a power of 2!\n");
    return;
  }

  if (vtxdist[npes] < (int)((float)(npes*npes)*1.2)) {
    if (mype == 0)
      printf("Error: Too many processors for this many vertices.\n");
    return;
  }

  minnvtxs = vtxdist[1]-vtxdist[0];
  for (i=0; i<npes; i++)
    minnvtxs = (minnvtxs < vtxdist[i+1]-vtxdist[i]) ? minnvtxs : vtxdist[i+1]-vtxdist[i];

  if (minnvtxs < (int)((float)npes*1.1)) {
    if (mype == 0)
      printf("Error: vertices are not distributed equally.\n");
    return;
  }
 
  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, order, npes, mype, 1);


  /*****************************/
  /* Set up control structures */
  /*****************************/
  if (options == NULL && options[0] == 0) {
    dbglvl_original = GLOBAL_DBGLVL;
    seed            = GLOBAL_SEED;
  }
  else {
    dbglvl_original = options[PMV3_OPTION_DBGLVL];
    seed            = options[PMV3_OPTION_SEED];
  }

  SetUpCtrl(&ctrl, nparts, 0, *comm);

  ctrl.CoarsenTo   = amin(vtxdist[npes]+1, 25*amax(npes, nparts));
  ctrl.seed        = (seed == 0 ? mype : seed*mype);
  ctrl.sync        = GlobalSEMax(&ctrl, seed);
  ctrl.partType    = STATIC_PARTITION;
  ctrl.ps_relation = -1;
  ctrl.tpwgts      = fsmalloc(nparts, 1.0/(float)(nparts), "tpwgts");
  ctrl.ubvec[0]    = 1.03;

  graph = Mc_SetUpGraph(&ctrl, 1, vtxdist, xadj, NULL, adjncy, NULL, &wgtflag);

  AllocateWSpace(&ctrl, graph, &wspace);

  /*=======================================================
   * Compute the initial k-way partitioning 
   =======================================================*/
  IFSET(dbglvl_original, DBG_TIME, InitTimers(&ctrl));
  IFSET(dbglvl_original, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(dbglvl_original, DBG_TIME, starttimer(ctrl.TotalTmr));

  Mc_Global_Partition(&ctrl, graph, &wspace);

  /* Collapse the number of partitions to be from 0..npes-1 */
  for (i=0; i<graph->nvtxs; i++)
    graph->where[i] = graph->where[i]%npes;
  ctrl.nparts = nparts = npes;

  /*=======================================================
   * Move the graph according to the partitioning
   =======================================================*/
  IFSET(dbglvl_original, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(dbglvl_original, DBG_TIME, starttimer(ctrl.MoveTmr));

  MALLOC_CHECK(NULL);
  graph->ncon = 1;
  mgraph = Mc_MoveGraph(&ctrl, graph, &wspace);
  MALLOC_CHECK(NULL);

  IFSET(dbglvl_original, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(dbglvl_original, DBG_TIME, stoptimer(ctrl.MoveTmr));

  /* restore the user supplied dbglvl */
  ctrl.dbglvl = dbglvl_original;

  /*=======================================================
   * Now compute an ordering of the moved graph
   =======================================================*/
  AdjustWSpace(&ctrl, mgraph, &wspace);

  ctrl.ipart = ISEP_NODE;
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, amax(20*npes, 1000));

  /* compute tvwgts */
  for (j=0; j<mgraph->ncon; j++)
    ltvwgts[j] = 0;

  for (i=0; i<mgraph->nvtxs; i++)
    for (j=0; j<mgraph->ncon; j++)
      ltvwgts[j] += mgraph->vwgt[i*mgraph->ncon+j];

  for (j=0; j<mgraph->ncon; j++)
    ctrl.tvwgts[j] = GlobalSESum(&ctrl, ltvwgts[j]);

  mgraph->nvwgt = fmalloc(mgraph->nvtxs*mgraph->ncon, "mgraph->nvwgt");
  for (i=0; i<mgraph->nvtxs; i++)
    for (j=0; j<mgraph->ncon; j++)
      mgraph->nvwgt[i*mgraph->ncon+j] = (float)(mgraph->vwgt[i*mgraph->ncon+j]) / (float)(ctrl.tvwgts[j]);


  morder = idxmalloc(mgraph->nvtxs, "PAROMETIS: morder");
  MultilevelOrder(&ctrl, mgraph, morder, sizes, &wspace);

  MALLOC_CHECK(NULL);

  /* Invert the ordering back to the original graph */
  ProjectInfoBack(&ctrl, graph, order, morder, &wspace);

  MALLOC_CHECK(NULL);

  IFSET(dbglvl_original, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(dbglvl_original, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(dbglvl_original, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(dbglvl_original, DBG_TIME, MPI_Barrier(ctrl.gcomm));

  GKfree((void **)&ctrl.tpwgts, &morder, LTERM);
  FreeGraph(mgraph);
  FreeInitialGraphAndRemap(graph, 0, 1);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, order, npes, mype, 0);

  MALLOC_CHECK(NULL);
}


/***********************************************************************************
* This function is the entry point of the parallel ordering algorithm.
* This function assumes that the graph is already nice partitioned among the 
* processors and then proceeds to perform recursive bisection.
************************************************************************************/
void PAROMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
                idxtype *order, idxtype *sizes, int *options, MPI_Comm comm)
{
  int numflag, newoptions[5];

  newoptions[0] = 1;
  newoptions[PMV3_OPTION_DBGLVL] = options[4];
  newoptions[PMV3_OPTION_SEED] = GLOBAL_SEED;

  numflag = options[3];

  ParMETIS_V3_NodeND(vtxdist, xadj, adjncy, &numflag, newoptions, order, sizes, &comm);

  options[0] = -1;

}

