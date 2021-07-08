C========================================================================
C     Broken up fortran interface to mumps solver, callable from C++.
C     Based on demo driver code in 
C
C       MUMPS 4.8.4, built on Mon Dec 15 15:31:38 UTC 2008
C
C========================================================================


c======================================================================
c     Module as common block replacement
c======================================================================
      module mumps_module

      INCLUDE 'dmumps_struc.h'

      INTEGER n_mumps
      INTEGER, dimension(:), allocatable :: workspace_scaling_factor
      TYPE (DMUMPS_STRUC), dimension(:), allocatable :: mumps_par_pool

      end module mumps_module


c======================================================================
c     Setup pool of mumps solvers
c======================================================================
      subroutine mumps_setup_solver_pool(max_n_solvers)

      use mumps_module

      INTEGER max_n_solvers

      n_mumps=max_n_solvers
      ALLOCATE( mumps_par_pool(max_n_solvers) )
      ALLOCATE( workspace_scaling_factor(max_n_solvers) )

      end


c======================================================================
c     Catch wrong ID
c======================================================================
      subroutine mumps_catch_wrong_id(i_pool)

      use mumps_module

      INTEGER i_pool

      if (i_pool .gt. n_mumps) then
         write(*,*) "i_pool ", i_pool," exceeds max. number of"
         write(*,*) "mumps solvers in pool (n_mumps=", n_mumps,")"
         stop
      endif

      end

c======================================================================
c     Setup routine
c======================================================================
      subroutine mumps_setup(i_pool,initial_workspace_scaling_factor)

      use mumps_module

      integer i_pool

      INTEGER IERR, I, nproc, nrows_per_proc
      integer initial_workspace_scaling_factor
      double precision error, max_error

c.....Range checking       
      call mumps_catch_wrong_id(i_pool)

C.....Initialize an instance of the package
C.....for L U factorization (sym = 0, with working host)

C.....Initialise     
      mumps_par_pool(i_pool)%JOB = -1
 
C.....Matrix is not symmetric
      mumps_par_pool(i_pool)%SYM = 0

C.....We're computing on the host too
      mumps_par_pool(i_pool)%PAR = 1

C.....Do it...      
      CALL DMUMPS(mumps_par_pool(i_pool))

C.....Output stream for global info on host. Negative value suppresses printing
      mumps_par_pool(i_pool)%icntl(3)=-1
      
C.....Only show error messages and stats
      mumps_par_pool(i_pool)%icntl(4)=2

C.....Assembled matrix (rather than element-by_element)
      mumps_par_pool(i_pool)%icntl(5)=0
      
C.....Distributed problem with user-specified distribution
      mumps_par_pool(i_pool)%icntl(18)=3
      
C.....Dense RHS
      mumps_par_pool(i_pool)%icntl(20)=0

C.....Non-distributed solution
      mumps_par_pool(i_pool)%icntl(21)=0

C.....Default scaling factor for workspace
      workspace_scaling_factor(i_pool)=initial_workspace_scaling_factor

      return
      end




c======================================================================
c     Specify multiplier for workspace scaling factor. Default is 2.
c======================================================================
      subroutine mumps_set_workspace_scaling_factor(i_pool,s)

      use mumps_module 

      integer i_pool
      double precision s

c.....Range checking       
      call mumps_catch_wrong_id(i_pool)

C.....Assign...
      workspace_scaling_factor(i_pool)=s

      return
      end



c======================================================================
c     Switch on doc of stats etc. Redirects output to "channel 6"
c======================================================================
      subroutine mumps_switch_on_doc(i_pool)

      use mumps_module

      integer i_pool

c.....Range checking       
      call mumps_catch_wrong_id(i_pool)

C.....Output stream for global info on host. Negative value suppresses printing
      mumps_par_pool(i_pool)%icntl(3)=6

      return
      end



c======================================================================
c     Switch off doc of stats etc. 
c======================================================================
      subroutine mumps_switch_off_doc(i_pool)


      use mumps_module

      integer i_pool

c.....Range checking       
      call mumps_catch_wrong_id(i_pool)

C.....Output stream for global info on host. Negative value suppresses printing
      mumps_par_pool(i_pool)%icntl(3)=-1

      return
      end

     


c======================================================================
c     Solver routine
c======================================================================
      subroutine mumps_solve(i_pool,n,nz_loc,irn_loc,jcn_loc,a_loc,rhs)

      use mumps_module

      integer i_pool

      integer n, nz_loc
      integer irn_loc(nz_loc), jcn_loc(nz_loc)
      integer initial_workspace_scaling_factor
      double precision a_loc(nz_loc), rhs(n)

c.....Range checking       
      call mumps_catch_wrong_id(i_pool)
      
C.....call mumps_setup
      initial_workspace_scaling_factor=2
      call mumps_setup(i_pool,initial_workspace_scaling_factor)
      
C.....analyse/factorise
      call mumps_factorise(i_pool,n,nz_loc,irn_loc,jcn_loc,a_loc)
      
C.....back-substitute
      call mumps_backsub(i_pool,n,rhs)

C.....Shutdown      
      CALL mumps_shutdown(i_pool)

      return
      end


c======================================================================
c     Cleanup memory routine
c======================================================================
      subroutine mumps_cleanup_memory(i_pool)

      use mumps_module

      integer i_pool

c.....Range checking       
      call mumps_catch_wrong_id(i_pool)

C.....Deallocate user data
      IF ( associated(mumps_par_pool(i_pool)%IRN_loc) )THEN
         DEALLOCATE( mumps_par_pool(i_pool)%IRN_loc)
      end if
      IF ( associated(mumps_par_pool(i_pool)%JCN_loc) )THEN         
         DEALLOCATE( mumps_par_pool(i_pool)%JCN_loc)
      end if
      IF ( associated(mumps_par_pool(i_pool)%A_loc) )THEN                  
         DEALLOCATE( mumps_par_pool(i_pool)%A_loc)
      end if
      IF ( associated(mumps_par_pool(i_pool)%RHS) )THEN                           
         DEALLOCATE( mumps_par_pool(i_pool)%RHS)
      END IF
  
      RETURN

      END




c======================================================================
c     Shutdown mumps
c======================================================================
      subroutine mumps_shutdown(i_pool)

      use mumps_module

      integer i_pool

c.....Range checking       
      call mumps_catch_wrong_id(i_pool)

C.....Cleanup memory just to be on the safe side
      CALL mumps_cleanup_memory(i_pool)

C.....Destroy the instance (deallocate internal data structures)
      mumps_par_pool(i_pool)%JOB = -2
      CALL DMUMPS(mumps_par_pool(i_pool))
  
      RETURN

      END


c======================================================================
c     Analyse and factorise
c======================================================================
      subroutine mumps_factorise(i_pool,n,nz_loc,irn_loc,jcn_loc,a_loc)

      use mumps_module

      integer i_pool

      integer n, nz_loc
      integer irn_loc(nz_loc),jcn_loc(nz_loc)
      double precision a_loc(nz_loc)

      INTEGER IERR, I, nproc, nrows_per_proc
      double precision error, max_error

c.....Range checking       
      call mumps_catch_wrong_id(i_pool)

C.....Specify size of system
      mumps_par_pool(i_pool)%N=n
      mumps_par_pool(i_pool)%NZ_loc=nz_loc

C.....Allocate storage for entries
      ALLOCATE( mumps_par_pool(i_pool)%IRN_loc( 
     &                      mumps_par_pool(i_pool)%NZ_loc ) )
      ALLOCATE( mumps_par_pool(i_pool)%JCN_loc( 
     &                      mumps_par_pool(i_pool)%NZ_loc ) )
      ALLOCATE( mumps_par_pool(i_pool)%A_loc( 
     &                      mumps_par_pool(i_pool)%NZ_loc ) )

      DO I = 1, mumps_par_pool(i_pool)%NZ_loc
         mumps_par_pool(i_pool)%IRN_loc(I)=irn_loc(I)
         mumps_par_pool(i_pool)%JCN_loc(I)=jcn_loc(I)
         mumps_par_pool(i_pool)%A_loc(I)=a_loc(I)
      END DO


C.....Do analysis
      mumps_par_pool(i_pool)%JOB = 1
      CALL DMUMPS(mumps_par_pool(i_pool))

C.....Document estimate for size of workspace
      if (mumps_par_pool(i_pool)%icntl(3).ge.0) then
         if (mumps_par_pool(i_pool)%myid.eq.0) then
            write(*,*) 'Estimated max. workspace in MB: ', 
     &           mumps_par_pool(i_pool)%infog(26)
         end if
      end if

c.....Jump address for re-factorisation with larger workspace
 10   continue

C.....Set workspace to multiple of that -- ought to be "significantly
C.....larger than infog(26)", according to the manual :(
      mumps_par_pool(i_pool)%icntl(23)=
     &    workspace_scaling_factor(i_pool)*
     &    mumps_par_pool(i_pool)%infog(26)

C.....Do factorisation
      mumps_par_pool(i_pool)%JOB = 2
      CALL DMUMPS(mumps_par_pool(i_pool))
  
C.....check for error
      if (mumps_par_pool(i_pool)%infog(1).ne.0) then
C........Increase scaling factor for workspace and run again         
         if (mumps_par_pool(i_pool)%myid.eq.0) then
            write(*,*) 'Processor', mumps_par_pool(i_pool)%myid, 
     &                 'has error during mumps factorisation!'
            write(*,*) 'Processor', mumps_par_pool(i_pool)%myid, 
     &                 'error codes: ',
     &                  mumps_par_pool(i_pool)%info(1), " ",
     &                  mumps_par_pool(i_pool)%info(2)
         end if
         workspace_scaling_factor(i_pool)=
     &   workspace_scaling_factor(i_pool)*2
         if (mumps_par_pool(i_pool)%myid.eq.0) then
            write(*,*) 'Increasing workspace_scaling_factor to ', 
     &                  workspace_scaling_factor(i_pool)
         end if
         go to 10
      else
         if (mumps_par_pool(i_pool)%icntl(3).ge.0) then
            if (mumps_par_pool(i_pool)%myid.eq.0) then
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
      subroutine mumps_backsub(i_pool,n,rhs)

      use mumps_module

      integer i_pool

      INTEGER IERR, I, nproc, nrows_per_proc
      double precision rhs(n)


c.....Range checking       
      call mumps_catch_wrong_id(i_pool)

      ALLOCATE( mumps_par_pool(i_pool)%RHS( n ) )
      DO I=1,n
         mumps_par_pool(i_pool)%RHS(I)=rhs(I)
      END DO

C.....Do solve/back-substitution
      mumps_par_pool(i_pool)%JOB = 3
      CALL DMUMPS(mumps_par_pool(i_pool))

C.....copy it back 
      DO I=1,n
         rhs(I)=mumps_par_pool(i_pool)%RHS(I)
      END DO
  
      RETURN

      END

