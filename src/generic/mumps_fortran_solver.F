C========================================================================
C     Broken up fortran interface to mumps solver, callable from C++.
C     Based on demo driver code in 
C
C       MUMPS 4.8.4, built on Mon Dec 15 15:31:38 UTC 2008
C
C========================================================================



c======================================================================
c     Setup routine
c======================================================================
      subroutine mumps_setup()

      INCLUDE 'dmumps_struc.h'

      INTEGER workspace_scaling_factor
      TYPE (DMUMPS_STRUC) mumps_par
      common /stuff/mumps_par,workspace_scaling_factor

      INTEGER IERR, I, nproc, nrows_per_proc
      double precision error, max_error
      
C.....Initialize an instance of the package
C.....for L U factorization (sym = 0, with working host)

C.....Initialise     
      mumps_par%JOB = -1
 
C.....Matrix is not symmetric
      mumps_par%SYM = 0

C.....We're computing on the host too
      mumps_par%PAR = 1

C.....Do it...      
      CALL DMUMPS(mumps_par)

C.....Output stream for global info on host. Negative value suppresses printing
      mumps_par%icntl(3)=-1
      
C.....Only show error messages and stats
      mumps_par%icntl(4)=2

C.....Assembled matrix (rather than element-by_element)
      mumps_par%icntl(5)=0
      
C.....Distributed problem with user-specified distribution
      mumps_par%icntl(18)=3
      
C.....Dense RHS
      mumps_par%icntl(20)=0

C.....Non-distributed solution
      mumps_par%icntl(21)=0

C.....Default scaling factor for workspace
      workspace_scaling_factor=1

      return
      end




c======================================================================
c     Specify multiplier for workspace scaling factor. Default is 2.
c======================================================================
      subroutine mumps_set_workspace_scaling_factor(s)

      INCLUDE 'dmumps_struc.h'

      INTEGER workspace_scaling_factor
      TYPE (DMUMPS_STRUC) mumps_par
      common /stuff/mumps_par,workspace_scaling_factor

      double precision s

C.....Assign...
      workspace_scaling_factor=s

      return
      end



c======================================================================
c     Switch on doc of stats etc. Redirects output to "channel 6"
c======================================================================
      subroutine mumps_switch_on_doc()

      INCLUDE 'dmumps_struc.h'

      INTEGER workspace_scaling_factor
      TYPE (DMUMPS_STRUC) mumps_par
      common /stuff/mumps_par,workspace_scaling_factor


C.....Output stream for global info on host. Negative value suppresses printing
      mumps_par%icntl(3)=6

      return
      end



c======================================================================
c     Switch off doc of stats etc. 
c======================================================================
      subroutine mumps_switch_off_doc()

      INCLUDE 'dmumps_struc.h'

      INTEGER workspace_scaling_factor
      TYPE (DMUMPS_STRUC) mumps_par
      common /stuff/mumps_par,workspace_scaling_factor


C.....Output stream for global info on host. Negative value suppresses printing
      mumps_par%icntl(3)=-1

      return
      end

     


c======================================================================
c     Solver routine
c======================================================================
      subroutine mumps_solve(n,nz_loc,irn_loc,jcn_loc,a_loc,rhs)

      INCLUDE 'dmumps_struc.h'

      INTEGER workspace_scaling_factor
      TYPE (DMUMPS_STRUC) mumps_par
      common /stuff/mumps_par,workspace_scaling_factor

      integer n, nz_loc
      integer irn_loc(nz_loc), jcn_loc(nz_loc)
      double precision a_loc(nz_loc), rhs(n)
      
C.....call mumps_setup
      call mumps_setup()
      
C.....analyse/factorise
      call mumps_factorise(n,nz_loc,irn_loc,jcn_loc,a_loc)
      
C.....back-substitute
      call mumps_backsub(n,rhs)

C.....Shutdown      
      CALL mumps_shutdown()

      return
      end


c======================================================================
c     Cleanup memory routine
c======================================================================
      subroutine mumps_cleanup_memory()

      INCLUDE 'dmumps_struc.h'

      INTEGER workspace_scaling_factor
      TYPE (DMUMPS_STRUC) mumps_par
      common /stuff/mumps_par,workspace_scaling_factor

C.....Deallocate user data
      IF ( associated(mumps_par%IRN_loc) )THEN
         DEALLOCATE( mumps_par%IRN_loc)
      end if
      IF ( associated(mumps_par%JCN_loc) )THEN         
         DEALLOCATE( mumps_par%JCN_loc)
      end if
      IF ( associated(mumps_par%A_loc) )THEN                  
         DEALLOCATE( mumps_par%A_loc)
      end if
      IF ( associated(mumps_par%RHS) )THEN                           
         DEALLOCATE( mumps_par%RHS)
      END IF
  
      RETURN

      END




c======================================================================
c     Shutdown mumps
c======================================================================
      subroutine mumps_shutdown()

      INCLUDE 'dmumps_struc.h'

      INTEGER workspace_scaling_factor
      TYPE (DMUMPS_STRUC) mumps_par
      common /stuff/mumps_par,workspace_scaling_factor

C.....Cleanup memory just to be on the safe side
      CALL mumps_cleanup_memory()

C.....Destroy the instance (deallocate internal data structures)
      mumps_par%JOB = -2
      CALL DMUMPS(mumps_par)
  
      RETURN

      END


c======================================================================
c     Analyse and factorise
c======================================================================
      subroutine mumps_factorise(n,nz_loc,irn_loc,jcn_loc,a_loc)

      INCLUDE 'dmumps_struc.h'

      INTEGER workspace_scaling_factor
      TYPE (DMUMPS_STRUC) mumps_par
      common /stuff/mumps_par,workspace_scaling_factor

      integer n, nz_loc
      integer irn_loc(nz_loc),jcn_loc(nz_loc)
      double precision a_loc(nz_loc)

      INTEGER IERR, I, nproc, nrows_per_proc
      double precision error, max_error

C.....Specify size of system
      mumps_par%N=n
      mumps_par%NZ_loc=nz_loc


C.....Allocate storage for entries
      ALLOCATE( mumps_par%IRN_loc( mumps_par%NZ_loc ) )
      ALLOCATE( mumps_par%JCN_loc( mumps_par%NZ_loc ) )
      ALLOCATE( mumps_par%A_loc( mumps_par%NZ_loc ) )

      DO I = 1, mumps_par%NZ_loc
         mumps_par%IRN_loc(I)=irn_loc(I)
         mumps_par%JCN_loc(I)=jcn_loc(I)
         mumps_par%A_loc(I)=a_loc(I)
      END DO


C.....Do analysis
      mumps_par%JOB = 1
      CALL DMUMPS(mumps_par)

C.....Document estimate for size of workspace
      if (mumps_par%icntl(3).ge.0) then
         if (mumps_par%myid.eq.0) then
            write(*,*) 'Estimated max. workspace in MB: ', 
     &           mumps_par%infog(26)
         end if
      end if

c.....Jump address for re-factorisation with larger workspace
 10   continue

C.....Set workspace to multiple of that -- ought to be "significantly
C.....larger than infog(26)", according to the manual :(
      mumps_par%icntl(23)=workspace_scaling_factor*mumps_par%infog(26)

C.....Do factorisation
      mumps_par%JOB = 2
      CALL DMUMPS(mumps_par)
  
C.....check for error
      if (mumps_par%infog(1).ne.0) then
C........Increase scaling factor for workspace and run again         
         if (mumps_par%myid.eq.0) then
            write(*,*) 'Error during mumps factorisation!'
         end if
         workspace_scaling_factor=workspace_scaling_factor*2
         if (mumps_par%myid.eq.0) then
            write(*,*) 'Increasing workspace_scaling_factor to ', 
     &                  workspace_scaling_factor
         end if
         go to 10
      else
         if (mumps_par%icntl(3).ge.0) then
            if (mumps_par%myid.eq.0) then
               write(*,*) 
     &         'Successfully completed factorisation in mumps'
            end if
         end if
      end if


      RETURN

      END


c======================================================================
c     Back-subsititute/solve
c======================================================================
      subroutine mumps_backsub(n,rhs)

      INCLUDE 'dmumps_struc.h'

      INTEGER workspace_scaling_factor
      TYPE (DMUMPS_STRUC) mumps_par
      common /stuff/mumps_par,workspace_scaling_factor

      INTEGER IERR, I, nproc, nrows_per_proc
      double precision rhs(n)

      ALLOCATE( mumps_par%RHS( n ) )
      DO I=1,n
         mumps_par%RHS(I)=rhs(I)
      END DO


C.....Do solve/back-substitution
      mumps_par%JOB = 3
      CALL DMUMPS(mumps_par)

C.....copy it back 
      DO I=1,n
         rhs(I)=mumps_par%RHS(I)
      END DO
  
      RETURN

      END

