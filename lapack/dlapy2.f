! Copyright (c) 1992-2013 The University of Tennessee and The University
!                         of Tennessee Research Foundation.  All rights
!                         reserved.
! Copyright (c) 2000-2013 The University of California Berkeley. All
!                         rights reserved.
! Copyright (c) 2006-2013 The University of Colorado Denver.  All rights
!                         reserved.
!
! $COPYRIGHT$
!
! Additional copyrights may follow
!
! $HEADER$
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
! - Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
!
! - Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer listed
!   in this license in the documentation and/or other materials
!   provided with the distribution.
!
! - Neither the name of the copyright holders nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
! The copyright holders provide no reassurances that the source code
! provided does not infringe any patent, copyright, or any other
! intellectual property rights of third parties.  The copyright holders
! disclaim any liability to any recipient for claims brought against
! recipient by any third party for infringement of that parties
! intellectual property rights.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*> \brief \b DLAPY2 returns sqrt(x2+y2).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLAPY2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   X, Y
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*> overflow.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is DOUBLE PRECISION
*>          X and Y specify the values x and y.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2017
*
*> \ingroup OTHERauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
      LOGICAL            X_IS_NAN, Y_IS_NAN
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN
      EXTERNAL           DISNAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      X_IS_NAN = DISNAN( X )
      Y_IS_NAN = DISNAN( Y )
      IF ( X_IS_NAN ) DLAPY2 = X
      IF ( Y_IS_NAN ) DLAPY2 = Y
*
      IF ( .NOT.( X_IS_NAN.OR.Y_IS_NAN ) ) THEN
         XABS = ABS( X )
         YABS = ABS( Y )
         W = MAX( XABS, YABS )
         Z = MIN( XABS, YABS )
         IF( Z.EQ.ZERO ) THEN
            DLAPY2 = W
         ELSE
            DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
         END IF
      END IF
      RETURN
*
*     End of DLAPY2
*
      END
