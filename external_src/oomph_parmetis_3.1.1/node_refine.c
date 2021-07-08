/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * node_refine.c
 *
 * This file contains code that performs the k-way refinement
 *
 * Started 3/1/96
 * George
 *
 * $Id: node_refine.c,v 1.2 2003/07/21 17:18:50 karypis Exp $
 */

#include <parmetislib.h>


/************************************************************************************/
/*! 
  This function allocates the memory required for the nodeND refinement code.
  The only refinement-related information that is has is the \c graph->where vector
  and allocates the memory for the remaining of the refinement-related data-structures.

*/
/************************************************************************************/
void AllocateNodePartitionParams(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int nparts, nvtxs;
  idxtype *vwgt;
  NRInfoType *rinfo, *myrinfo;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayInitTmr));

  nvtxs  = graph->nvtxs;
  nparts = ctrl->nparts;

  graph->nrinfo  = (NRInfoType *)GKmalloc(sizeof(NRInfoType)*nvtxs, "AllocateNodePartitionParams: rinfo");
  graph->lpwgts  = idxmalloc(2*nparts, "AllocateNodePartitionParams: lpwgts");
  graph->gpwgts  = idxmalloc(2*nparts, "AllocateNodePartitionParams: gpwgts");
  graph->sepind  = idxmalloc(nvtxs, "AllocateNodePartitionParams: sepind");
  graph->hmarker = idxmalloc(nvtxs, "AllocateNodePartitionParams: hmarker");

  /* Allocate additional memory for graph->vwgt in order to store the weights
     of the remote vertices */
  vwgt        = graph->vwgt;
  graph->vwgt = idxmalloc(nvtxs+graph->nrecv, "AllocateNodePartitionParams: graph->vwgt");
  idxcopy(nvtxs, vwgt, graph->vwgt);
  GKfree((void **)&vwgt, LTERM);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayInitTmr));
}


/************************************************************************************/
/*! 
  This function computes the initial node refinment information for the parallel
  nodeND code. It requires that the required data-structures have already been
  allocated via a call to AllocateNodePartitionParams.

*/
/************************************************************************************/
void ComputeNodePartitionParams(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, j, nparts, nvtxs, nsep, firstvtx, lastvtx;
  idxtype *xadj, *adjncy, *adjwgt, *vtxdist, *vwgt, *lpwgts, *gpwgts, *sepind, *hmarker;
  idxtype *where;
  NRInfoType *rinfo, *myrinfo;
  int me, other, otherwgt;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayInitTmr));

  nvtxs  = graph->nvtxs;
  nparts = ctrl->nparts;

  vtxdist = graph->vtxdist;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  vwgt    = graph->vwgt;

  where   = graph->where;
  rinfo   = graph->nrinfo;
  lpwgts  = graph->lpwgts;
  gpwgts  = graph->gpwgts;
  sepind  = graph->sepind;
  hmarker = graph->hmarker;

  firstvtx = vtxdist[ctrl->mype];
  lastvtx  = vtxdist[ctrl->mype+1];

  /* Reset refinement data structures */
  idxset(2*nparts, 0, lpwgts);
  idxset(nvtxs, 0, hmarker);

  /* Send/Receive the where and vwgt information of interface vertices. */
  CommInterfaceData(ctrl, graph, where, wspace->indices, where+nvtxs); 
  CommInterfaceData(ctrl, graph, vwgt, wspace->indices, vwgt+nvtxs); 


  /*------------------------------------------------------------
  / Compute now the degrees
  /------------------------------------------------------------*/
  for (nsep=i=0; i<nvtxs; i++) {
    me = where[i];
    ASSERT(ctrl, me >= 0 && me < 2*nparts);
    lpwgts[me] += vwgt[i];

    if (me >= nparts) {  /* If it is a separator vertex */
      sepind[nsep++] = i;
      lpwgts[2*nparts-1] += vwgt[i];  /* Keep track of total separator weight */

      myrinfo = rinfo+i;
      myrinfo->edegrees[0] = myrinfo->edegrees[1] = 0;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (me != other)
          myrinfo->edegrees[other%2] += vwgt[adjncy[j]];
      }
    }
  }
  graph->nsep = nsep;

  /* Finally, sum-up the partition weights */
  MPI_Allreduce((void *)lpwgts, (void *)gpwgts, 2*nparts, IDX_DATATYPE, MPI_SUM, ctrl->comm);
  graph->mincut = gpwgts[2*nparts-1];


  /* Mark the halo vertices by determining how many non-local vertices they are
   * connected to. */
  idxset(nvtxs, 0, hmarker);
  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (adjncy[j] >= nvtxs) {
        hmarker[i]++;
      }
    }
  }


#ifdef XX
  /* Print Weight information */
  if (ctrl->mype == 0) {
    for (i=0; i<nparts; i+=2) 
      printf("[%5d %5d %5d] ", gpwgts[i], gpwgts[i+1], gpwgts[nparts+i]); 
    printf("\n");
  }
#endif

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayInitTmr));
}


/************************************************************************************/
/*! 
  This function performs k-way node-based refinement. The refinement is done
  concurrently for all the different partitions. It works because each of the
  partitions is disconnected from each other due to the removal of the previous level
  separators.

*/
/************************************************************************************/
void KWayNodeRefine(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace, 
         int npasses, float ubfraction)
{
  int i, ii, iii, j, k, pass, nvtxs, firstvtx, lastvtx, otherlastvtx, c, nmoves, 
      nlupd, nsupd, nnbrs, nchanged, nsep;
  int npes = ctrl->npes, mype = ctrl->mype, nparts = ctrl->nparts;
  idxtype *xadj, *adjncy, *adjwgt, *vtxdist, *vwgt;
  idxtype *where, *lpwgts, *gpwgts, *sepind;
  idxtype *peind, *recvptr, *sendptr;
  idxtype *update, *supdate, *rupdate, *pe_updates, *htable, *changed;
  idxtype *badminpwgt, *badmaxpwgt;
  KeyValueType *swchanges, *rwchanges;
  int *nupds_pe;
  NRInfoType *rinfo, *myrinfo;
  int from, to, me, other, oldcut;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayTmr));

  nvtxs = graph->nvtxs;

  vtxdist = graph->vtxdist;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  vwgt    = graph->vwgt;

  firstvtx = vtxdist[mype];
  lastvtx  = vtxdist[mype+1];

  where  = graph->where;
  rinfo  = graph->nrinfo;
  lpwgts = graph->lpwgts;
  gpwgts = graph->gpwgts;

  nsep   = graph->nsep;
  sepind = graph->sepind;

  nnbrs   = graph->nnbrs;
  peind   = graph->peind;
  recvptr = graph->recvptr;
  sendptr = graph->sendptr;

  changed   = idxmalloc(nvtxs, "KWayRefine: changed");
  rwchanges = wspace->pairs;
  swchanges = rwchanges + recvptr[nnbrs];

  update   = idxmalloc(nvtxs, "KWayRefine: update");
  supdate  = wspace->indices;
  rupdate  = supdate + recvptr[nnbrs];
  nupds_pe = imalloc(npes, "KWayRefine: nupds_pe");

  htable = idxsmalloc(nvtxs+graph->nrecv, 0, "KWayRefine: lhtable");

  badminpwgt = wspace->pv1;
  badmaxpwgt = wspace->pv2;

  for (i=0; i<nparts; i+=2) {
    badminpwgt[i] = badminpwgt[i+1] = (1.0/ubfraction)*(gpwgts[i]+gpwgts[i+1])/2;
    badmaxpwgt[i] = badmaxpwgt[i+1] = ubfraction*(gpwgts[i]+gpwgts[i+1])/2;
  }
  //myprintf(ctrl, "%6d %6d %6d %6d %6d %6d %6d\n", lpwgts[0], lpwgts[1], lpwgts[2], gpwgts[0], gpwgts[1], gpwgts[2], badmaxpwgt[0]);

  IFSET(ctrl->dbglvl, DBG_REFINEINFO, 
      PrintNodeBalanceInfo(ctrl, nparts, gpwgts, badminpwgt, badmaxpwgt, 1));

  for (pass=0; pass<npasses; pass++) {
    oldcut = graph->mincut;

    for (c=0; c<2; c++) {
      for (i=0; i<nparts; i+=2) {
        badminpwgt[i] = badminpwgt[i+1] = (1.0/ubfraction)*(gpwgts[i]+gpwgts[i+1])/2;
        badmaxpwgt[i] = badmaxpwgt[i+1] = ubfraction*(gpwgts[i]+gpwgts[i+1])/2;
      }

      nlupd = nsupd = nmoves = nchanged = 0;
      for (ii=0; ii<nsep; ii++) {
        i = sepind[ii];
        from = where[i];

        ASSERT(ctrl, from >= nparts);

        /* Go through the loop to see if gain is possible for the separator vertex */
        if (rinfo[i].edegrees[(c+1)%2] <= vwgt[i]) {
          /* It is a one-sded move so it will go to the other partition. 
             Look at the comments in InitMultisection to understand the meaning 
             of from%nparts */
          to = from%nparts+c;  

          if (gpwgts[to]+vwgt[i] > badmaxpwgt[to]) {
            /* printf("Skip because of weight! %d\n", vwgt[i]-rinfo[i].edegrees[(c+1)%2]); */
            continue;   /* We cannot move it there because it gets too heavy */
          }

          /* Update the where information of the vertex you moved */
          where[i] = to;

          /* Remove this vertex from the sepind. Note the trick for looking at 
             the sepind[ii] again */
          sepind[ii--] = sepind[--nsep]; 

          /* myprintf(ctrl, "Vertex %d [%d %d] is moving to %d from %d [%d]\n", 
                  i+firstvtx, vwgt[i], rinfo[i].edegrees[(c+1)%2], to, from, where[i]); */

          lpwgts[from]       -= vwgt[i];
          lpwgts[2*nparts-1] -= vwgt[i];
          lpwgts[to]         += vwgt[i];
          gpwgts[to]         += vwgt[i];

          /* Put the vertices adjacent to i that belong to either the separator or
             the (c+1)%2 partition into the update array */
          for (j=xadj[i]; j<xadj[i+1]; j++) {
            k = adjncy[j];
            if (htable[k] == 0 && where[k] != to) {
              htable[k] = 1;
              if (k<nvtxs)
                update[nlupd++] = k;
              else
                supdate[nsupd++] = k;
            }
          }
          nmoves++;
          if (graph->pexadj[i+1]-graph->pexadj[i] > 0)
            changed[nchanged++] = i;
        }
      }

      /* myprintf(ctrl, "nmoves: %d, nlupd: %d, nsupd: %d\n", nmoves, nlupd, nsupd); */

      /* Tell everybody interested what the new where[] info is for the interface vertices */
      CommChangedInterfaceData(ctrl, graph, nchanged, changed, where, swchanges, 
          rwchanges, wspace->pv4); 


      IFSET(ctrl->dbglvl, DBG_RMOVEINFO, rprintf(ctrl, "\t[%d %d], [%d %d %d]\n", 
                pass, c, GlobalSESum(ctrl, nmoves), GlobalSESum(ctrl, nsupd), 
                GlobalSESum(ctrl, nlupd)));


      /*-----------------------------------------------------------------------
      / Time to communicate with processors to send the vertices whose degrees 
      / need to be updated.
      /-----------------------------------------------------------------------*/
      /* Issue the receives first */
      for (i=0; i<nnbrs; i++) {
        MPI_Irecv((void *)(rupdate+sendptr[i]), sendptr[i+1]-sendptr[i], IDX_DATATYPE,
                  peind[i], 1, ctrl->comm, ctrl->rreq+i);
      }

      /* Issue the sends next. This needs some preporcessing */
      for (i=0; i<nsupd; i++) {
        htable[supdate[i]] = 0;
        supdate[i] = graph->imap[supdate[i]];
      }
      iidxsort(nsupd, supdate);

      for (j=i=0; i<nnbrs; i++) {
        otherlastvtx = vtxdist[peind[i]+1];
        for (k=j; k<nsupd && supdate[k] < otherlastvtx; k++); 
        MPI_Isend((void *)(supdate+j), k-j, IDX_DATATYPE, peind[i], 1, ctrl->comm, 
            ctrl->sreq+i);
        j = k;
      }

      /* OK, now get into the loop waiting for the send/recv operations to finish */
      MPI_Waitall(nnbrs, ctrl->rreq, ctrl->statuses);
      for (i=0; i<nnbrs; i++) 
        MPI_Get_count(ctrl->statuses+i, IDX_DATATYPE, nupds_pe+i);
      MPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses);


      /*-------------------------------------------------------------
      / Place the received to-be updated vertices into update[] 
      /-------------------------------------------------------------*/
      for (i=0; i<nnbrs; i++) {
        pe_updates = rupdate+sendptr[i];
        for (j=0; j<nupds_pe[i]; j++) {
          k = pe_updates[j];
          if (htable[k-firstvtx] == 0) {
            htable[k-firstvtx] = 1;
            update[nlupd++] = k-firstvtx;
          }
        }
      }


      /*-------------------------------------------------------------
      / Update the where information of the vertices that are pulled
      / into the separator.
      /-------------------------------------------------------------*/
      nchanged = 0;
      for (ii=0; ii<nlupd; ii++) {
        i = update[ii];
        me = where[i];
        if (me < nparts && me%2 == (c+1)%2) { /* This vertex is pulled into the separator */
          lpwgts[me] -= vwgt[i];
          where[i] = nparts+me-(me%2); 
          sepind[nsep++] = i;  /* Put the vertex into the sepind array */
          if (graph->pexadj[i+1]-graph->pexadj[i] > 0)
            changed[nchanged++] = i;

          lpwgts[where[i]]   += vwgt[i];
          lpwgts[2*nparts-1] += vwgt[i];
          /* myprintf(ctrl, "Vertex %d moves into the separator from %d to %d\n", 
                 i+firstvtx, me, where[i]); */
        }
      }

      /* Tell everybody interested what the new where[] info is for the interface vertices */
      CommChangedInterfaceData(ctrl, graph, nchanged, changed, where, swchanges, 
          rwchanges, wspace->pv4); 


      /*-------------------------------------------------------------
      / Update the rinfo of the vertices in the update[] array
      /-------------------------------------------------------------*/
      for (ii=0; ii<nlupd; ii++) {
        i = update[ii];
        ASSERT(ctrl, htable[i] == 1);

        htable[i] = 0;

        me = where[i];
        if (me >= nparts) {  /* If it is a separator vertex */
          /* myprintf(ctrl, "Updating %d %d\n", i+firstvtx, me); */

          myrinfo = rinfo+i;
          myrinfo->edegrees[0] = myrinfo->edegrees[1] = 0;

          for (j=xadj[i]; j<xadj[i+1]; j++) {
            other = where[adjncy[j]];
            if (me != other)
              myrinfo->edegrees[other%2] += vwgt[adjncy[j]];
          }
        }
      }

      /* Finally, sum-up the partition weights */
      MPI_Allreduce((void *)lpwgts, (void *)gpwgts, 2*nparts, IDX_DATATYPE, MPI_SUM, 
          ctrl->comm);
      graph->mincut = gpwgts[2*nparts-1];

      IFSET(ctrl->dbglvl, DBG_REFINEINFO, PrintNodeBalanceInfo(ctrl, nparts, gpwgts, 
            badminpwgt, badmaxpwgt, 0));
    }

    if (graph->mincut == oldcut)
      break;
  }

  GKfree((void **)&update, &nupds_pe, &htable, &changed, LTERM);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayTmr));
}



/*************************************************************************
* This function prints balance information for the parallel k-section 
* refinement algorithm
**************************************************************************/
void PrintNodeBalanceInfo(CtrlType *ctrl, int nparts, idxtype *gpwgts, idxtype *badminpwgt, 
         idxtype *badmaxpwgt, int title)
{
  int i;

  if (ctrl->mype == 0) {
    if (title)
      printf("K-way sep-refinement: TotalSep: %d, ", gpwgts[2*nparts-1]);
    else
      printf("\tTotalSep: %d, ", gpwgts[2*nparts-1]);

    for (i=0; i<nparts; i+=2) 
      printf(" [%5d %5d %5d %5d %5d]", gpwgts[i], gpwgts[i+1], gpwgts[nparts+i], 
          badminpwgt[i], badmaxpwgt[i]); 
    printf("\n");
  }
  MPI_Barrier(ctrl->comm);
}

