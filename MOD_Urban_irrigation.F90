#include <define.h>

MODULE MOD_Urban_Irrigation

  ! -----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! 城市植被中添加灌溉部分
  ! -----------------------------------------------------------------------
  ! !USE

USE MOD_Precision
USE MOD_Vars_Global
USE MOD_TimeManager, only: julian2monthday, isleapyear
USE MOD_Urban_LUCY, only: gmt2local
IMPLICIT NONE
SAVE


PUBLIC :: MOD_Urban_Irrigation



CONTAINS
!-----------------------------------------------------------------------


   SUBROUTINE Urban_Irrigation (nl_soil,        idate,         fveg         porsl, &
                                bsw,            psi0,          deltim,      dz_soisno, &
                                t_soisno,       wliq_soisno,   patchlonr)

  !=======================================================================
  ! 计算总的灌溉水量和当前时间步长的灌溉水量
  !=======================================================================

    use precision
    use MOD_Const_Physical, only : tfrz
    implicit none

  !-----------------------Argument-----------------------------------------

    Integer,  Intent(in) :: nl_soil                ! upper bound of array
    INTEGER,  Intent(in) :: idate(3)               ! calendar (year, julian day, seconds)
    real(r8), INTENT(in) :: fveg                   ! fraction of vegetation cover
    real(r8), INTENT(in) :: porsl(1:nl_soil)       ! soil porosity [-]
    real(r8), INTENT(in) :: bsw(1:nl_soil)         ! Clapp-Hornberger "B"
    real(r8), INTENT(in) :: psi0(1:nl_soil)        ! saturated soil suction (mm) (NEGATIVE)
    real(r8), INTENT(in) :: deltim                 ! seconds in a time step [second]
    real(r8), INTENT(in) :: dz_soisno(1:nl_soil)   ! layer thickness (m)
    real(r8), INTENT(in) :: t_soisno(1:nl_soil)    ! soil/snow skin temperature (K)
    real(r8), INTENT(in) :: wliq_soisno(1:nl_soil) ! liquid water (kg/m2)
    real(r8), INTENT(in) :: patchlonr              ! longitude of patch [radian]

    real(r8), INTENT(out) :: irr_layer(1:nl_soil)  ! irrigation water of a layer, all layers add to 1
    real(r8), INTENT(out) :: irr_total             ! total irrigation water
    real(r8), INTENT(out) :: irr_current           ! 当前时间步长的灌溉量 (mm h2o/s)
    

    real(r8)              :: londeg                ! longitude of path [degree]

  !-----------------------Local Variables------------------------------

    real(r8) s_node            ! vol_liq/porosity
    real(r8) smpmax            ! wilting point potential in mm
    real(r8) smpfc             ! Field capacity potential in mm
    real(r8) smp_node          ! matrix potential
    real(r8) t_total           ! 总的灌溉时间
    real(r8) s_target          ! 目标含水量

    INTEGER :: &
         ldate(3), &! local time (year, julian day, seconds)
         ihour   , &! hour of day
         day     , &! day of mmonth
         month   , &! month of year
         t_hour  , &! 灌溉时间
         i          ! loop counter

  !-----------------------End Variables list---------------------------

        ! Total amount of irrigation (irr_total) and 当前灌溉水量 (irr_current)

         
    if (fveg > 1.e-5 .and. LAI >1.e-5) then

        londeg = patchlonr*180/PI

        CALL gmt2local(idate, londeg, ldate)
        CALL julian2monthday(ldate(1), ldate(2), month, day)

        if (month==7 .or. month==8 .or. month==9) then

            ihour = CEILING(ldate(3)*1./3600)

                if (t_hour == ihour) then 
                    irr_total = 1.e-10 
                        do i = 1, nl_soil
                          if(t_soisno(i)>tfrz .and. porsl(i)>=1.e-6)then

                             smpmax     = -1.5e5
                             smpfc      = -3.4e3

                             s_node     = max(wliq_soisno(i)/(1000.*dz_soisno(i)),0.001)
                             s_node     = min(1., s_node)

                             s_target   = max(porsl(i)*(smpfc(i)/psi0(i)**(-bsw(i)),0.001) 
                             s_target   = min(1.,s_target)

                                 if (s_node     <=   s_target) then
                                    irr_layer(i) =   dz_soisno(i)*(s_target - s_node)
                                 endif  

                          else
                             irr_layer(i) = 0
                          endif

                        irr_total    =   irr_total + irr_layer(i)


                        end do
                    a_irr               =   irr_total/t_total
                    irr_current         =   a_irr *deltim  
                endif  
            endif
        endif
        
   end subroutine urban_irrigation


END MOD_Urban_Irrigation