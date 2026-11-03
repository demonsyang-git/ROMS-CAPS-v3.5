#include "cppdefs.h"
      MODULE ice_frazil_mod
#if (defined ICE_MODEL && defined ICE_THERMO && !defined TS_FIXED) \
      || defined CICE_MODEL
!
!=======================================================================
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the frazil ice growth in the water when the
!  water temperature gets below freezing. It adjusts the water
!  temperature and salinity accordingly.
!
!  Reference: Steele et al. (1989). JPO, 19, 139-147.
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC ice_frazil

      CONTAINS

      SUBROUTINE ice_frazil (ng, tile)

      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_ice
      USE mod_stepping
# ifdef CICE_MODEL
      USE mod_forces
# endif

      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 80, __LINE__, __FILE__)
# endif
!
      CALL ice_frazil_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nnew(ng),                                   &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % z_r,                             &
     &                      OCEAN(ng) % rho,                            &
     &                      OCEAN(ng) % t,                              &
# ifdef CICE_MODEL
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr,                         &
     &                      FORCES(ng) % stflx,                         &
     &                      ICE(ng) % ai,                               &
     &                      ICE(ng) % hi,                               &
     &                      GRID(ng) % z_w,                             &
# endif
     &                      ICE(ng) % wfr)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 80, __LINE__, __FILE__)
# endif
      RETURN
      END SUBROUTINE ice_frazil
!
!***********************************************************************
      subroutine ice_frazil_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nnew,                                 &
# ifdef MASKING
     &                            rmask,                                &
# endif
# ifdef WET_DRY
     &                            rmask_wet,                            &
# endif
     &                            Hz, z_r, rho, t,                      &
# ifdef CICE_MODEL
     &                            sustr, svstr, stflx, ai, hi, z_w,     &
# endif
     &                            wfr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# ifdef CICE_MODEL
      USE shr_const_mod
      use icepack_intfc, only: icepack_sea_freezing_temperature,        &
     &                         icepack_query_parameters
# endif
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod
      USE distribute_mod, ONLY : mp_reduce
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nnew

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#  ifdef CICE_MODEL
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
      real(r8), intent(in) :: ai(LBi:,LBj:)
      real(r8), intent(in) :: hi(LBi:,LBj:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  endif
      real(r8), intent(out) :: wfr(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  ifdef CICE_MODEL
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(out) :: wfr(LBi:UBi,LBj:UBj)
# endif
!
! Local variable definitions
!
      integer :: i, j, k, itrc

      real(r8), parameter :: Lhat = 79.2_r8
      real(r8), parameter :: r = 0.5_r8

      !real(r8) :: t_freeze !! cyy masked
      real(r8) :: s1
      real(r8) :: z1
      real(r8) :: sice
      real(r8) :: gamma_k
      real(r8) :: t_fr
      real(r8) :: ice_dens
      real(r8) :: delta_wfr
      real(r8), parameter :: depressT = -0.054_r8 !! same as CICE
      !real(r8), parameter :: depressT = -0.0543_r8

# ifdef CICE_MODEL
      real(r8) :: dz, heat, tfr
      real(r8), parameter :: bl_thick = 5.0_r8
      real(r8), parameter :: sl_thick = 1.0_r8
      real(r8) :: ice_thick, thickness, fract
      real(r8), parameter :: latent_heat_fusion = 3.337e5_r8   ! J/kg
      real(r8), parameter :: cp_sw              = 3996.0_r8    ! J/kg/K
      real(r8), parameter :: rho_sw             = 1026.0_r8    ! kg/m^3
      real(r8), parameter :: rho_fw             = 1000.0_r8    ! kg/m^3
      real(r8), parameter :: cp_over_lhfusion   =                       &
     &                 rho_sw*cp_sw/(latent_heat_fusion*rho_fw)
      real(r8), parameter :: rhocpr = 4093740._r8       ! [W s m-3 K-1]
      real(r8), parameter :: salice = 4.0_r8
      real(r8), parameter :: salref = 34.7_r8
      real(r8) :: potice, salocn !! cyy added
      logical :: update_ocn_f    !! cyy added
# endif

# ifdef DISTRIBUTE
      real(r8), allocatable           :: buffer(:)
      character (len=3), allocatable  :: op_handle(:)
# endif

!  Inline functions for freezing points, cyy note: using CICE function now
!  Freezing temperature (Gill, 1982)
!     t_freeze(s1,z1) = -0.0575*s1 + 1.710523d-3*sqrt(s1)**3
!    &       - 2.154996d-4*s1*s1 + 0.000753*z1
!  Freezing temperature (Steele et al. 1989)
!      t_freeze(s1,z1) = -0.0543*s1 + 0.000759*z1
!      t_freeze(s1,z1) = depressT*s1
!      t_freeze(s1,z1) = -1.8_r8

# include "set_bounds.h"

# ifdef DISTRIBUTE
      IF (.not. allocated(buffer))    allocate (buffer(1))
      IF (.not. allocated(op_handle)) allocate (op_handle(1))
# endif

      DO j=Jstr,Jend
        DO i=Istr,Iend
          wfr(i,j) = 0.0_r8
        ENDDO
      ENDDO

# ifdef CICE_MODEL
      call icepack_query_parameters(update_ocn_f_out=update_ocn_f)

      DO j=Jstr,Jend
        DO i=Istr,Iend
#  if defined WET_DRY
          IF (rmask_wet(i,j) .ne. 0.0_r8) THEN
#  endif
#  if defined MASKING
          IF (rmask(i,j) .ne. 0.0_r8) THEN
#  endif
            DO k=N(ng),1,-1
              dz = Hz(i,j,k)
              t_fr = icepack_sea_freezing_temperature(t(i,j,k,nnew,isalt))
              potice = (t_fr - t(i,j,k,nnew,itemp)) * dz

              !*** if potice < 0, compute again later
              !*** if potice > 0, keep on freezing (wfr > 0)
              IF ((wfr(i,j) + potice) > 0) THEN
                wfr(i,j) = wfr(i,j) + potice
              ELSE
                potice = -wfr(i,j)
                wfr(i,j) = 0.0_r8
              ENDIF
              t(i,j,k,nnew,itemp) = t(i,j,k,nnew,itemp)                 &
     &              + potice/Hz(i,j,k)
              if (.not.update_ocn_f) then
                t(i,j,k,nnew,isalt) = t(i,j,k,nnew,isalt)               &
     &              + (salref-salice)*potice*cp_over_lhfusion/Hz(i,j,k)
              endif
              !*** jump out if done with surface layer
              IF (z_w(i,j,N(ng))-z_w(i,j,k-1) > bl_thick) EXIT
            END DO
#  if defined MASKING
          END IF
#  endif
#  if defined WET_DRY
          END IF
#  endif
        END DO
      END DO

! Now compute heat flux from ocean and subtract from wfr.
      DO j=Jstr,Jend
        DO i=Istr,Iend
#  if defined WET_DRY
          IF (rmask_wet(i,j) .ne. 0.0_r8) THEN
#  endif
#  if defined MASKING
          IF (rmask(i,j) .ne. 0.0_r8) THEN
#  endif
            heat = 0.0_r8
            thickness = 0.0_r8
            DO k=N(ng),1,-1
              dz = Hz(i,j,k)
              tfr = icepack_sea_freezing_temperature(t(i,j,k,nnew,isalt))
              IF ((thickness + dz) < sl_thick) THEN
                thickness = thickness + dz
                heat = heat + (t(i,j,k,nnew,itemp)-tfr)*dz
              ELSE
                fract = (sl_thick - thickness)/dz
                !thickness = thickness + dz*fract
                heat = heat + (t(i,j,k,nnew,itemp)-tfr)*dz*fract
              END IF
              IF (z_w(i,j,N(ng))-z_w(i,j,k-1) > sl_thick) EXIT
            END DO
            ! convert units to W/m^2
            wfr(i,j) = (wfr(i,j) - heat)*rhocpr/dt(ng)
#  if defined MASKING
          END IF
#  endif
#  if defined WET_DRY
          END IF
#  endif
        END DO
      END DO
# endif

      CALL bc_r2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          wfr)
# ifdef DISTRIBUTE
      buffer(1) = exit_flag
      op_handle(1) = 'MAX'
      CALL mp_reduce (ng, iNLM, 1, buffer, op_handle)
      exit_flag = int(buffer(1))
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    wfr)
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,nnew,:))
# endif
      RETURN
      END SUBROUTINE ice_frazil_tile

#endif
      END MODULE ice_frazil_mod

