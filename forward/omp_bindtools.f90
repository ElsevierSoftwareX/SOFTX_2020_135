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

!>    @brief init thread binding table for ScaleMP systems
      SUBROUTINE init_scalemp_binding()
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     node-size: core number to use per node, offset for the next node
        INTEGER node_size, node_offs
        INTEGER o, i, bs, k, ni
        ! INTEGER g_o, g_i
!     binding table
        INTEGER po, pi, mo, mi
        PARAMETER (po=128)
        PARAMETER (pi=128)
        INTEGER (kind=4) btable(pi,po)
        COMMON /binding_tab/mo, mi, btable
!     node list
        INTEGER inode_list(po+1), nnode
        character (len=256) :: cnode_list
!     error handler
        INTEGER error
!
        INTRINSIC trim, mod
#ifdef libnuma
        INTEGER (kind=4) f_numa_num_configured_nodes
        INTEGER (kind=4) f_numa_num_configured_cpus
        INTEGER (kind=4) f_numa_available
        EXTERNAL f_numa_num_configured_nodes
        EXTERNAL f_numa_num_configured_cpus
        EXTERNAL f_numa_available
#endif

!       default for old ScaleMP system (on RWTH Aachen university)
        node_size = 8
        node_offs = node_size
#ifdef libnuma
!       overwrite for modern systems using libnuma
        IF (f_numa_available()==-1) THEN
          WRITE(*,'(1A)') &
            'error: "libnuma" not working, see "omp_bindtools.f90"!'
          STOP
        END IF
        node_size = f_numa_num_configured_cpus() / &
          f_numa_num_configured_nodes()
        node_offs = node_size
        WRITE(*,'(1A,1I4)') '  [I] : ccNUMA nodes = ', &
          f_numa_num_configured_nodes()
        WRITE(*,'(1A,1I4)') &
          '  [I] : cores per ccNUMA node = ', node_size
#endif
!
        CALL get_environment_variable(name='JOB_NODE_LIST', &
          value=cnode_list,status=error)
        IF (error==0) THEN
!        default value to overwrite
          DO i = 1, po + 1
            inode_list(i) = -1
          END DO
          WRITE(*,*) ' [I] : JOB_NODE_LIST = "', trim(cnode_list), &
            '"'
!        read node list
          READ(cnode_list,*,err=101,end=101) (inode_list(i),i=1,po)
101       nnode = 0
!        count node number
          DO i = 1, po
            IF (inode_list(i)>=0) nnode = nnode + 1
          END DO
        ELSE
!        default setup, assume all nodes
          WRITE(*,'(2A)') '  <D> : JOB_NODE_LIST not defined,', &
            ' using all ccNUMA nodes'
          nnode = 1
#ifdef libnuma
          nnode = f_numa_num_configured_nodes()
#endif
          DO i = 1, nnode
            inode_list(i) = i-1
          END DO
          inode_list(nnode+1) = -1
        END IF
        WRITE(*,'(1A,999I6)') '        nodes = ', &
          (inode_list(i),i=1,nnode)
!
        mi = tlevel_1
        mo = tlevel_0
!      threading block size [bs] for inner parallel region
        bs = (node_size*nnode) / tlevel_0
!
        IF (bs>=tlevel_1) THEN
          WRITE(*,'(1A,3I6)') '  [I] : binding mode = memory distance'
!        compute binding list
          DO o = 1, mo
            DO i = 1, mi
              k = (o-1)*bs +(i-1)*bs/tlevel_1
              ni = k / node_size + 1
              btable(i,o) = inode_list(ni)*node_offs +k-(ni-1)*node_size
            END DO
          END DO
        ELSE
          WRITE(*,*) ' [E] : not enough cores for all threads'
          STOP
        END IF
!     show table
        WRITE(*,*) ' [I] : binding table'
        WRITE(*,'(1A,2I6)') '        ', mo, mi
        DO o = 1, mo
          WRITE(*,'(1A,128I6)') '        ', (btable(i,o),i=1,mi)
        END DO
!
        RETURN
      END

!>    @brief load thread binding table
!>    @param[in] fname file name of the table
      SUBROUTINE load_binding(fname)
        IMPLICIT NONE
        character (len=*) :: fname
        INTEGER po, pi, mo, mi
#ifdef TBIND
        INTEGER o, i
#endif
        PARAMETER (po=128)
        PARAMETER (pi=128)
        INTEGER (kind=4) :: btable(pi,po)
        COMMON /binding_tab/mo, mi, btable
        INTRINSIC trim

#ifdef TBIND
        IF (fname=='default') THEN
          WRITE(*,*) ' <D> : no thread binding'
          mo = 0
          mi = 0
        ELSE
          WRITE(*,*) ' [R] : thread binding table from "', &
            trim(fname), '"'
          OPEN(81,file=fname,status='OLD',err=100)
          READ(81,*) mo, mi
          IF (mo>po .OR. mi>pi .OR. mo<=0 .OR. mi<=0) THEN
            WRITE(*,'(2(1A,1I4),1A)') &
              'error: processor binding table exceed limits (', po, &
              'x', pi, ') in "load_binding" !'
            STOP
          END IF
          DO o = 1, mo
            READ(81,*) (btable(i,o),i=1,mi)
          END DO
          CLOSE(81)
        END IF
!     show table
        WRITE(*,*) ' [I] : binding table'
        WRITE(*,'(1A,2I6)') '        ', mo, mi
        DO o = 1, mo
          WRITE(*,'(1A,128I6)') '        ', (btable(i,o),i=1,mi)
        END DO
#endif
        RETURN
#ifdef TBIND
100     WRITE(*,'(3A)') 'error: can not open map file "', &
          trim(fname), '" !'
#endif
        STOP
      END

!>    @brief try thread binding
!>    @param[in] o outer thread index
      SUBROUTINE omp_binding(o)
        use mod_genrl
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INTEGER o
#ifdef TBIND
        INTEGER lo
#endif
        INCLUDE 'OMP_TOOLS.inc'
        INTEGER p_o, p_i, mo, mi
        PARAMETER (p_o=128)
        PARAMETER (p_i=128)
        INTEGER (kind=4) :: btable(p_i,p_o)
        COMMON /binding_tab/mo, mi, btable
!     master binding staff
#ifdef TBIND
        INTEGER (kind=4) :: p_is
#endif
#ifdef libnuma
        INTEGER (kind=4) :: f_numa_run_on_node
        INTEGER (kind=4 :: f_numa_node_of_cpu
        INTEGER (kind=4) :: f_numa_available
        EXTERNAL f_numa_run_on_node
        EXTERNAL f_numa_node_of_cpu
        EXTERNAL f_numa_available
#endif
        INTRINSIC mod

#ifdef libnuma
        IF (f_numa_available()==-1) THEN
          WRITE(*,'(1A)') &
            'error: "libnuma" not working, see "omp_bindtools.f90"!'
          STOP
        END IF
#endif
#ifdef TBIND
!       binding enabled?
        IF (mo==0 .OR. mi==0) RETURN
!
        i = omp_get_his_thread_num() + 1
        lo = o
!       special handling for ENKF/SIMUL, when fewer threads than realisations
        IF (lo>Tlevel_0 .AND. lo<=nsmpl) lo = mod(lo-1,Tlevel_0) +1
!
        IF (lo>mo .OR. i>mi .OR. lo<1 .OR. i<1) THEN
          WRITE(*,'(2(1A,1I4),1A)') &
            'error: processor binding table exceed limits (', mo, 'x', &
            mi, ') in "omp_binding" !'
          WRITE(*,'(2(1A,1I4),1A)') &
            '       for processor binding on (', lo, 'x', i, ')'
          STOP
        END IF
#ifdef libnuma
        IF (f_numa_run_on_node(f_numa_node_of_cpu(btable(i,lo)))<0) THEN
          WRITE(*,'(1A)') &
            'error: "libnuma" not working, see "omp_bindtools.f90"!'
          STOP
        END IF
#else
#ifndef G95
!     Intel, PGI, SUN
        CALL r_processorbind(btable(i,lo))
#else
!     GNU 4.xx, to avoid THREAD-spinning on the master affinity
        IF (i>1) THEN
          CALL r_processorbind(btable(i,lo))
        ELSE
!        (GNU compiler workaround) bind the master thread on all team cores (group binding)
          p_is = mi
          CALL r_groupbind(p_is,btable(1,lo))
        END IF
#endif
#endif
#endif
        RETURN
      END
