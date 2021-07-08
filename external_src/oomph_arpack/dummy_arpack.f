      SUBROUTINE BREAKAR()
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) '===================================================='
      WRITE(*,*) ' '
      WRITE(*,*) 'Broken dummy code! Please obtain the full'
      WRITE(*,*) 'ARPACK Fortran sources from '
      WRITE(*,*) ' '
      WRITE(*,*) '   http://www.caam.rice.edu/software/ARPACK/ '
      WRITE(*,*) ' '
      WRITE(*,*) 'Place all the ARPACK routines in '
      WRITE(*,*) ' '
      WRITE(*,*) '    external_src/oomph_arpack/all_arpack_sources.f'
      WRITE(*,*) ' '
      WRITE(*,*) 'and re-run '
      WRITE(*,*) ' '
      WRITE(*,*) '     autogen.sh --rebuild'
      WRITE(*,*) ' '
      WRITE(*,*) 'in the top-level oomph-lib directory.'
      WRITE(*,*) ' '
      WRITE(*,*) '[Terminating code execution] '
      WRITE(*,*) ' '
      WRITE(*,*) '==================================================='
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      STOP
      END

C*******************************************************************

c HERE ARE THE BROKEN STUBS FOR THE SUBROUTINES REQUIRED FROM HSL


C*******************************************************************
      subroutine dnaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, 
     &     ipntr, workd, workl, lworkl, info )

      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol

      integer    iparam(11), ipntr(14)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)

      call BREAKAR()
      end

      subroutine dneupd (rvec , howmny, select, dr    , di,    
     &                   z    , ldz   , sigmar, sigmai, workev,
     &                   bmat , n     , which , nev   , tol,
     &                   resid, ncv   , v     , ldv   , iparam,
     &                   ipntr, workd , workl , lworkl, info)

      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision      
     &           sigmar, sigmai, tol
c
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision 
     &           dr(nev+1)    , di(nev+1), resid(n)  , 
     &           v(ldv,ncv)   , z(ldz,*) , workd(3*n), 
     &           workl(lworkl), workev(3*ncv)

      call BREAKAR()
      end
