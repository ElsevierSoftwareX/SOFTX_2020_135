! MIT License
!
! Copyright (c) 2020 SHEMAT-Suite
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

!>    @brief check parameter / section - definition
!>    @param[in] filename model file name
      SUBROUTINE read_check(filename)
        use mod_genrl
        IMPLICIT NONE
        integer :: k
        character (len=80) :: filename
        character (len=80) :: line
        LOGICAL phi_style
        INTEGER ifil, locstr, f_entry, lblank, line_count
        EXTERNAL locstr, lblank
        INTEGER i_default
        INTEGER (kind=8) :: i_64, test
!     Integer compatibility test for BLAS and LAPACK
!       both must be from the same size as the default integer size
        i_default = 131072
        i_64 = 131072
        i_default = i_default*i_default
        i_64 = i_64*i_64
        IF (i_default==i_64) THEN
          test = 0
#ifdef BLAS32
          test = 1
#endif
#ifdef LAPACK32
          test = 1
#endif
          IF (test==1) THEN
            WRITE(*,'(2A)') 'Integer compatibility test for', &
              ' BLAS and LAPACK library: [failure]'
            WRITE(*,*) &
              'Please determine the current integer size type'
            WRITE(*,*) '  of your BLAS and LAPACK library, both must'
            WRITE(*,*) &
              '  have the same size and equal to the one used'
            WRITE(*,*) '  by the compiler for this program binary !'
            STOP
          END IF
        END IF
        WRITE(*,*)
        WRITE(*,*) '  checking input definitions:'
        WRITE(*,*) '    from file "', filename(:lblank(filename)), &
          '"'
!     init counter for hard errors (old style - not longer supported)
        k = 0
!     read file
        phi_style = .FALSE.
        ifil = 77
        line_count = 0
        OPEN(ifil,file=filename,status='old')
10      CONTINUE
        READ(ifil,'(A)',end=20) line
        line_count = line_count + 1
        IF (line(1:1)/='#') GO TO 10
        f_entry = 0
!     forward, general
        f_entry = f_entry + locstr(line,key_char//' title')
        f_entry = f_entry + locstr(line,key_char//' PROPS')
        f_entry = f_entry + locstr(line,key_char//' USER')
        f_entry = f_entry + locstr(line,key_char//' linfo')
        f_entry = f_entry + locstr(line,key_char//' runmode')
        f_entry = f_entry + locstr(line,key_char//' grid')
        f_entry = f_entry + locstr(line,key_char//' nlsolve')
        f_entry = f_entry + locstr(line,key_char//' active')
        f_entry = f_entry + locstr(line,key_char//' rref')
        f_entry = f_entry + locstr(line,key_char//' grav')
        f_entry = f_entry + locstr(line,key_char//' rhocm')
        f_entry = f_entry + locstr(line,key_char//' hpf')
        f_entry = f_entry + locstr(line,key_char//' delx')
        f_entry = f_entry + locstr(line,key_char//' dely')
        f_entry = f_entry + locstr(line,key_char//' delz')
        f_entry = f_entry + locstr(line,key_char//' uindex')
        f_entry = f_entry + locstr(line,key_char//' units')
        f_entry = f_entry + locstr(line,key_char//' bcunits')
        f_entry = f_entry + locstr(line,key_char//' split units')
!     forward, output
        f_entry = f_entry + locstr(line,key_char//' borehole log')
        f_entry = f_entry + locstr(line,key_char//' disable small output')
        f_entry = f_entry + locstr(line,key_char//' disable output')
        f_entry = f_entry + locstr(line,key_char//' write output')
        f_entry = f_entry + locstr(line,key_char//' read output')
!     forward, data
        f_entry = f_entry + locstr(line,key_char//' data')
!     solver, convergency, flow
        f_entry = f_entry + locstr(line,key_char//' lsolvef')
        f_entry = f_entry + locstr(line,key_char//' error lsolvef')
        f_entry = f_entry + locstr(line,key_char//' maxiter lsolvef')
        f_entry = f_entry + locstr(line,key_char//' name lsolvef')
        f_entry = f_entry + locstr(line,key_char//' criteria lsolvef')
        f_entry = f_entry + locstr(line,key_char//' precondition lsolvef')
        f_entry = f_entry + locstr(line,key_char//' omf')
        f_entry = f_entry + locstr(line,key_char//' nliterf')
        f_entry = f_entry + locstr(line,key_char//' grad nliterf')
!     solver, convergency, temperature
        f_entry = f_entry + locstr(line,key_char//' lsolvet')
        f_entry = f_entry + locstr(line,key_char//' error lsolvet')
        f_entry = f_entry + locstr(line,key_char//' maxiter lsolvet')
        f_entry = f_entry + locstr(line,key_char//' name lsolvet')
        f_entry = f_entry + locstr(line,key_char//' criteria lsolvet')
        f_entry = f_entry + locstr(line,key_char//' precondition lsolvet')
        f_entry = f_entry + locstr(line,key_char//' omt')
        f_entry = f_entry + locstr(line,key_char//' nlitert')
        f_entry = f_entry + locstr(line,key_char//' grad nlitert')
!     solver, convergency, concentration
        f_entry = f_entry + locstr(line,key_char//' lsolvec')
        f_entry = f_entry + locstr(line,key_char//' error lsolvec')
        f_entry = f_entry + locstr(line,key_char//' maxiter lsolvec')
        f_entry = f_entry + locstr(line,key_char//' name lsolvec')
        f_entry = f_entry + locstr(line,key_char//' criteria lsolvec')
        f_entry = f_entry + locstr(line,key_char//' precondition lsolvec')
        f_entry = f_entry + locstr(line,key_char//' nliterc')
        f_entry = f_entry + locstr(line,key_char//' grad nliterc')
!     debug special
#ifdef DEBUG
        f_entry = f_entry + locstr(line,key_char//' debug output times')
#endif
!
!     forward, state variables
        f_entry = f_entry + locstr(line,key_char//' head init')
        f_entry = f_entry + locstr(line,key_char//' temp init')
        f_entry = f_entry + locstr(line,key_char//' pres init')
!     forward, boundary conditions
        f_entry = f_entry + locstr(line,key_char//' head bcd')
        f_entry = f_entry + locstr(line,key_char//' head bcn')
        f_entry = f_entry + locstr(line,key_char//' pres bcd')
        f_entry = f_entry + locstr(line,key_char//' pres bcn')
        f_entry = f_entry + locstr(line,key_char//' temp bcd')
        f_entry = f_entry + locstr(line,key_char//' temp bcn')
        f_entry = f_entry + locstr(line,key_char//' conc bcd')
        f_entry = f_entry + locstr(line,key_char//' conc bcn')
!     forward, time
        f_entry = f_entry + locstr(line,key_char//' timestep control')
        f_entry = f_entry + locstr(line,key_char//' tunit')
        f_entry = f_entry + locstr(line,key_char//' tstart')
        f_entry = f_entry + locstr(line,key_char//' titer')
        f_entry = f_entry + locstr(line,key_char//' time periods')
        f_entry = f_entry + locstr(line,key_char//' bc time periods')
        f_entry = f_entry + locstr(line,key_char//' monitor')
        f_entry = f_entry + locstr(line,key_char//' output times')
        f_entry = f_entry + locstr(line,key_char//' file output')
!     transport (chemical)
        f_entry = f_entry + locstr(line,key_char//' ntrans')
        f_entry = f_entry + locstr(line,key_char//' tracer')
        f_entry = f_entry + locstr(line,key_char//' transpar')
!     simulation
        f_entry = f_entry + locstr(line,key_char//' simulate')
        f_entry = f_entry + locstr(line,key_char//' samples')
        f_entry = f_entry + locstr(line,key_char//' parameter group')
!     ENKF
        f_entry = f_entry + locstr(line,key_char//' enkf postcompute')
        f_entry = f_entry + locstr(line,key_char//' enkf iter')
!     AD deterministic inverse
        f_entry = f_entry + locstr(line,key_char//' inverse')
        f_entry = f_entry + locstr(line,key_char//' enable property')
        f_entry = f_entry + locstr(line,key_char//' enable unit')
!     inverse, property units
        f_entry = f_entry + locstr(line,key_char//' errors')
        f_entry = f_entry + locstr(line,key_char//' apriori')
        f_entry = f_entry + locstr(line,key_char//' weight property')
        f_entry = f_entry + locstr(line,key_char//' optimize property')
!     inverse, boundary conditions
        f_entry = f_entry + locstr(line,key_char//' bcerrors')
        f_entry = f_entry + locstr(line,key_char//' bcapriori')
        f_entry = f_entry + locstr(line,key_char//' weight bc')
        f_entry = f_entry + locstr(line,key_char//' optimize bc')
!     inverse, time
        f_entry = f_entry + locstr(line,key_char//' tperrors')
        f_entry = f_entry + locstr(line,key_char//' tpapriori')
        f_entry = f_entry + locstr(line,key_char//' optimize tp')
!     inverse, weighting
!AW      f_entry = f_entry + locstr(line,key_char//' weight tp')
        f_entry = f_entry + locstr(line,key_char//' weight para')
        f_entry = f_entry + locstr(line,key_char//' weight data')
! --------
!
!     automatic generated check sections, by Makfile 'check_section'
       f_entry = f_entry + locstr(line,key_char//' set')
        f_entry = f_entry + locstr(line,key_char//' velocity')
        f_entry = f_entry + locstr(line,key_char//' mfd')
        f_entry = f_entry + locstr(line,key_char//' h5parse data file')
        f_entry = f_entry + locstr(line,key_char//' simtime')
        f_entry = f_entry + locstr(line,key_char//' jacoby')
        f_entry = f_entry + locstr(line,key_char//' itimestep')
        f_entry = f_entry + locstr(line,key_char//' iter_inv')
        f_entry = f_entry + locstr(line,key_char//' variable step size')
        f_entry = f_entry + locstr(line,key_char//' nliters')
        f_entry = f_entry + locstr(line,key_char//' grad nliters')
        f_entry = f_entry + locstr(line,key_char//' cindex')
        f_entry = f_entry + locstr(line,key_char//' simul')
        f_entry = f_entry + locstr(line,key_char//' enkf')
        f_entry = f_entry + locstr(line,key_char//' prop limit')
        f_entry = f_entry + locstr(line,key_char//' ilu block size')
        f_entry = f_entry + locstr(line,key_char//' ZETA')
        f_entry = f_entry + locstr(line,key_char//' THETA')
        f_entry = f_entry + locstr(line,key_char//' PSI')
        f_entry = f_entry + locstr(line,key_char//' LAMDA')
        f_entry = f_entry + locstr(line,key_char//' C0')
        f_entry = f_entry + locstr(line,key_char//' B2')
        f_entry = f_entry + locstr(line,key_char//' B1')
        f_entry = f_entry + locstr(line,key_char//' B0')
        f_entry = f_entry + locstr(line,key_char//' species')
        f_entry = f_entry + locstr(line,key_char//' secondary')
        f_entry = f_entry + locstr(line,key_char//' sec_comp')
        f_entry = f_entry + locstr(line,key_char//' parameter species')
        f_entry = f_entry + locstr(line,key_char//' minerals')
        f_entry = f_entry + locstr(line,key_char//' kinetic rate laws')
        f_entry = f_entry + locstr(line,key_char//' inhibitors')
        f_entry = f_entry + locstr(line,key_char//' subsample')
        f_entry = f_entry + locstr(line,key_char//' bc standard deviation')
        f_entry = f_entry + locstr(line,key_char//' bc group')
        f_entry = f_entry + locstr(line,key_char//' regular')
        f_entry = f_entry + locstr(line,key_char//' standard deviation')
        f_entry = f_entry + locstr(line,key_char//' covar prior para')
        f_entry = f_entry + locstr(line,key_char//' covar prior data')
        f_entry = f_entry + locstr(line,key_char//' fluid props')
        f_entry = f_entry + locstr(line,key_char//' archie')
        f_entry = f_entry + locstr(line,key_char//' gas props')
        f_entry = f_entry + locstr(line,key_char//' rhod')
        f_entry = f_entry + locstr(line,key_char//' liq init')
        f_entry = f_entry + locstr(line,key_char//' freeze_b')
        f_entry = f_entry + locstr(line,key_char//' freeze_a')
        f_entry = f_entry + locstr(line,key_char//' k_freeze')
        f_entry = f_entry + locstr(line,key_char//' fracpar')
! --------
        IF (f_entry==0) WRITE(*,'(A,I6,3A)') 'warning: ignore line ', &
          line_count, ': "', line(:lblank(line)), '"'
!     compatibility check
        IF (locstr(line,key_char//' phi')==1) phi_style = .TRUE.
! !!! no longer supported section, stoping later !!!
        k = k + locstr(line,key_char//' left head')
        k = k + locstr(line,key_char//' left flow')
        k = k + locstr(line,key_char//' right head')
        k = k + locstr(line,key_char//' right flow')
        k = k + locstr(line,key_char//' front head')
        k = k + locstr(line,key_char//' front flow')
        k = k + locstr(line,key_char//' back head')
        k = k + locstr(line,key_char//' back flow')
        k = k + locstr(line,key_char//' base head')
        k = k + locstr(line,key_char//' base flow')
        k = k + locstr(line,key_char//' top head')
        k = k + locstr(line,key_char//' top flow')
        k = k + locstr(line,key_char//' left temperature')
        k = k + locstr(line,key_char//' left heat flow')
        k = k + locstr(line,key_char//' right temperature')
        k = k + locstr(line,key_char//' right heat flow')
        k = k + locstr(line,key_char//' front temperature')
        k = k + locstr(line,key_char//' front heat flow')
        k = k + locstr(line,key_char//' back temperature')
        k = k + locstr(line,key_char//' back heat flow')
        k = k + locstr(line,key_char//' base temperature')
        k = k + locstr(line,key_char//' base heat flow')
        k = k + locstr(line,key_char//' top temperature')
        k = k + locstr(line,key_char//' top heat flow')
        GO TO 10
!     close project config file
20      CLOSE(ifil)
        IF (k>0) THEN
          WRITE(*,'(1A)') &
            'error: old BC sections no longer supported!!!'
          STOP
        END IF
!     compatibility check fails, then message out
        IF (phi_style) THEN
          WRITE(*,*)
          WRITE(*,'(1A)') &
            'error: found old "PHI" style definitions !!!'
          WRITE(*,'(1A)') &
            '       Please change all "phi" into "head" occurrence !'
          WRITE(*,*)
          STOP
        END IF
        RETURN
      END
