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

!> @brief calculate temperature dependent thermal conductivity
!> @param[in] lammref thermal conductivity from input file
!> @param[in] tlocal temperature
!> @param[in] tref reference temperature
!> @param[in] ismpl local sample index
!> @return temperature dependent thermal conductivity
!> @details
!> calculate temperature dependent thermal conductivity\n
!>
!>        lam_zoth/haenel = (770degC/(350degC+T) + 0.7) W/mK
!>
!> If T > 800degC, use the formula from zoth & haenel, 1988. \n
!>
!>        lamm = lam_zoth/haenel
!>
!> If T < 800degC, use the same formula with a factor `fct`
!>
!>        lamm = fct * lam_zoth/haenel
!>
!> The factor `fct` introduces an additional temperature dependence
!> such that\n
!> 1. `lamm(20degC) = lammref`\n
!> 2. `lamm(800degC) = lam_zoth/haenel(800degC)\n\n
!>
!> Thus, the thermal conductivity in the input file should resemble
!> information about the value of the thermal conductivity of the
!> matrix at temperature 20 degC.\n\n
!>
!> Sources: \n
!> Zoth, G., & Haenel, R, Handbook of terrestrial heat-flow density
!> determination, (1988), Appendix 10.1 Thermal
!> Conductivity. http://dx.doi.org/10.1007/978-94-009-2847-3\n\n
!>
!> Lehmann, H., Wang, K., & Clauser, C., Parameter identification and
!> uncertainty analysis for heat transfer at the ktb drill site using
!> a 2-d inverse method, Tectonophysics, 291(1-4), 179–194
!> (1998). Section 2.2 http://dx.doi.org/10.1016/s0040-1951(98)00039-0
!> \n
      double precision function lamm(lammref,tlocal,tref,ismpl)

        implicit none

        ! Sample index
        integer :: ismpl

        ! Thermal conductivity of matrix at tref=20degC
        double precision, intent (in) :: lammref

        ! Local temperature [degC]
        double precision, intent (in) :: tlocal

        ! Reference temperature [20 degC]
        double precision, intent (in) :: tref

        ! Upper limit temperature, where the approximation becomes
        ! equal to the general Zoth/Haenel formula
        double precision, parameter :: tlimit = 800.0d0

        ! lam_zoth/haenel at tlocal
        double precision :: lamm_zh

        ! lam_zoth/haenel at tref
        double precision :: lamm_zhref

        ! Reference Interpolation-factor
        double precision :: fctref

        ! Weight: Quotient of temperature differences
        double precision :: twgt

        ! Interpolation factor
        double precision :: fct


        if (tlocal > tlimit) then

          ! lam_zoth/haenel at tlocal
          lamm = 770.0d0/(350.0d0+tlocal) + 0.7D0

        else

          ! lam_zoth/haenel at tlocal
          lamm_zh = 770.0d0/(350.0d0+tlocal) + 0.7d0

          ! lam_zoth/haenel at tref
          lamm_zhref = 770.0d0/(350.0d0+tref) + 0.7d0

          ! Reference Interpolation-factor: Input lamm at tref divided
          ! by lam_zoth/haenel at tref
          fctref = lammref/lamm_zhref

          ! Quotient of temperature differences, local minus reference
          ! divided by limit minus reference
          twgt = (tlocal-tref)/(tlimit-tref)

          ! Interpolation factor between fctref at tref and 1 at tlimit
          fct = fctref*(1-twgt) + twgt

          ! Final lam: input at tref and lam_zoth/haenel at tlimit
          lamm = fct*lamm_zh

        end if

        return

      end function lamm
