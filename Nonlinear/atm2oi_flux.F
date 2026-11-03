#include "cppdefs.h"
      MODULE atm2oi_flux_mod
#ifdef ATM2OI_FLUXES
!
!svn $Id: atm2oi_flux.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== John C. Warner   ===
!                                                   Hernan G. Arango   =
!                                                   Chao-Yuan Yang     =
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the flux over the ocean surface               !
!  assuming all cells are ice-free, ice-ocean flux will be             !
!  weighted in the total flux afterward                                !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: atm2oi_flux

      CONTAINS
!
!***********************************************************************
      SUBROUTINE atm2oi_flux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_ocean
      USE mod_ice       !! cyy for the ice-weighted fluxes
      USE mod_grid
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 17)
# endif
      CALL atm2oi_flux_tile (ng, tile,                                  &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs(ng),                                 &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef WET_DRY
     &                        GRID(ng) % rmask_wet,                     &
     &                        GRID(ng) % umask_wet,                     &
     &                        GRID(ng) % vmask_wet,                     &
# endif
# ifdef EMINUSP
     &                        FORCES(ng) % EminusP,                     &
     &                        FORCES(ng) % evap,                        &
# endif
     &                        FORCES(ng) % rain,                        &
# ifdef SNOWFALL
     &                        FORCES(ng) % snow,                        &
     &                        FORCES(ng) % rain_lq,                     &
# endif
# ifdef WW3_MODEL
                              FORCES(ng) % Tauicex,                     &
                              FORCES(ng) % Tauicey,                     &
# endif               
     &                        FORCES(ng) % lhflx,                       &
     &                        OCEAN(ng)  % t,                           &
     &                        FORCES(ng) % lrflx,                       &
     &                        FORCES(ng) % nlrflx,                      &
     &                        FORCES(ng) % shflx,                       &
     &                        FORCES(ng) % snhflx,                      &
     &                        FORCES(ng) % albedo,                      &
     &                        FORCES(ng) % srflx,                       &
     &                        FORCES(ng) % nsrflx,                      &
     &                        FORCES(ng) % stflx,                       &
     &                        FORCES(ng) % LW_down,                     &
     &                        FORCES(ng) % SW_down,                     &
     &                        FORCES(ng) % Zair,                        &
     &                        FORCES(ng) % PotT,                        &
     &                        FORCES(ng) % Tair,                        &
     &                        FORCES(ng) % Qair,                        &
     &                        FORCES(ng) % Uwind,                       &
     &                        FORCES(ng) % Vwind,                       &
     &                        FORCES(ng) % sustr,                       &
     &                        FORCES(ng) % svstr,                       &
     &                        FORCES(ng) % Taux,                        &
     &                        FORCES(ng) % Tauy)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17)
# endif
      RETURN
      END SUBROUTINE atm2oi_flux
!
!***********************************************************************
      SUBROUTINE atm2oi_flux_tile (ng, tile,                            &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nrhs,                               &
# ifdef MASKING
     &                              rmask, umask, vmask,                &
# endif
# ifdef WET_DRY
     &                              rmask_wet, umask_wet, vmask_wet,    &
# endif
# ifdef EMINUSP
     &                              EminusP, evap,                      &
# endif
     &                              rain,                               &
# ifdef SNOWFALL
     &                              snow, rain_lq,                      &
# endif
# ifdef WW3_MODEL
                                    Tauicex, Tauicey,                   &
# endif                      
     &                              lhflx, t, lrflx, nlrflx,            &
     &                              shflx, snhflx,                      &
     &                              albedo, srflx, nsrflx, stflx,       &
     &                              LW_down, SW_down,                   &
     &                              zair, pott, tair, qair,             &
     &                              uwind, vwind,                       &
     &                              sustr, svstr,                       &
     &                              Taux, Tauy)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#  ifdef EMINUSP
      real(r8), intent(out) :: EminusP(LBi:,LBj:)
      real(r8), intent(inout) :: evap(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: rain(LBi:,LBj:)
#  ifdef SNOWFALL
      real(r8), intent(inout) :: snow(LBi:,LBj:)
      real(r8), intent(out) :: rain_lq(LBi:,LBj:)  !! cyy
#  endif
#  ifdef WW3_MODEL
      real(r8), intent(inout) :: Tauicex(LBi:,LBj:)
      real(r8), intent(inout) :: Tauicey(LBi:,LBj:)
#  endif      
      real(r8), intent(inout) :: lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
      real(r8), intent(out)   :: nlrflx(LBi:,LBj:) !! cyy
      real(r8), intent(inout) :: shflx(LBi:,LBj:)
      real(r8), intent(inout) :: snhflx(LBi:,LBj:)
      real(r8), intent(in)    :: albedo(LBi:,LBj:) !! cyy
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(out)   :: nsrflx(LBi:,LBj:) !! cyy
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(out)   :: LW_down(LBi:,LBj:)
      real(r8), intent(out)   :: SW_down(LBi:,LBj:)
      real(r8), intent(inout) :: zair(LBi:,LBj:)
      real(r8), intent(inout) :: pott(LBi:,LBj:)
      real(r8), intent(inout) :: tair(LBi:,LBj:)
      real(r8), intent(inout) :: qair(LBi:,LBj:)
      real(r8), intent(inout) :: uwind(LBi:,LBj:)
      real(r8), intent(inout) :: vwind(LBi:,LBj:)
      real(r8), intent(inout) :: sustr(LBi:,LBj:)
      real(r8), intent(inout) :: svstr(LBi:,LBj:)
      real(r8), intent(in)    :: Taux(LBi:,LBj:)
      real(r8), intent(in)    :: Tauy(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  ifdef EMINUSP
      real(r8), intent(out) :: EminusP(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: evap(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: rain(LBi:UBi,LBj:UBj)
#  ifdef SNOWFALL
      real(r8), intent(inout) :: snow(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rain_lq(LBi:UBi,LBj:UBj)   !! cyy
#  endif
      real(r8), intent(inout) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(out)   :: nlrflx(LBi:UBi,LBj:UBj)  !! cyy
      real(r8), intent(inout) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: snhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: albedo(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(out)   :: nsrflx(LBi:UBi,LBj:UBj)  !! cyy
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(out)   :: LW_down(LBi:UBi,LBj:UBj)
      real(r8), intent(out)   :: SW_down(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: zair(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: pott(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tair(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: qair(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vwind(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Taux(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Tauy(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff
      real(r8) :: Hscale, Hscale2
      real(r8) :: TseaK, LRad

# include "set_bounds.h"
!
!=======================================================================
!  Atmosphere-Ocean flux computations.
!=======================================================================
!
      cff=1.0_r8/rhow
      Hscale=rho0*Cp
      Hscale2=1.0_r8/(rho0*Cp)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ! CICE required variables
          ! downward radiation, C m/s to W/m2
          LW_down(i,j) = lrflx(i,j)*Hscale
          SW_down(i,j) = srflx(i,j)*Hscale
# ifdef SNOWFALL
          ! seperate liquid precip. from total precip.
          ! WRF precip. are all in water equivalent units 
          rain_lq(i,j) = rain(i,j) - snow(i,j)
          rain_lq(i,j) = max(0.0_r8, rain_lq(i,j)) !! cyy: need check
# endif
          ! net longwave radiation
          ! lrflx, nlrflx in C m/s
          TseaK = t(i,j,N(ng),nrhs,itemp) + 273.16_r8
!          LRad = lrflx(i,j)*Hscale-                                     &
!     &           emmiss*StefBo*TseaK*TseaK*TseaK*TseaK
          ! cyy note: the absorption of downward longwave radition
          !           should also consider the sea surface emissivity 
          LRad = emmiss*(lrflx(i,j)*Hscale-                             &
     &                   StefBo*TseaK*TseaK*TseaK*TseaK)
          nlrflx(i,j) = LRad*Hscale2
          ! net shortwave radiation
          ! srflx, nsrflx in C m/s
          nsrflx(i,j) = srflx(i,j)*(1._r8-albedo(i,j))

          ! net heat flux over the ocean surface
          snhflx(i,j)=(nsrflx(i,j) + nlrflx(i,j) +                      &
     &                 lhflx(i,j)  + shflx(i,j))
# ifdef EMINUSP
          stflx(i,j,isalt)=cff*(evap(i,j)-rain(i,j))
# endif
# ifdef MASKING
          snhflx(i,j)=snhflx(i,j)*rmask(i,j)
          rain(i,j)=rain(i,j)*rmask(i,j)
#   ifdef EMINUSP
          evap(i,j)=evap(i,j)*rmask(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask(i,j)
#   endif
# endif
# ifdef WET_DRY
          snhflx(i,j)=snhflx(i,j)*rmask_wet(i,j)
          rain(i,j)=rain(i,j)*rmask_wet(i,j)
#   ifdef EMINUSP
          evap(i,j)=evap(i,j)*rmask_wet(i,j)
          stflx(i,j,isalt)=stflx(i,j,isalt)*rmask_wet(i,j)
#   endif
# endif

# ifdef EMINUSP
          EminusP(i,j)=stflx(i,j,isalt)
# endif
        END DO
      END DO
!
!  Transfer Taux and Tauy from rho points to u and v points.
!
      cff=0.5_r8
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=cff*(Taux(i-1,j)+Taux(i,j))
# ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
# endif
# ifdef WET_DRY
          sustr(i,j)=sustr(i,j)*umask_wet(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=cff*(Tauy(i,j-1)+Tauy(i,j))
# ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
# endif
# ifdef WET_DRY
          svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
# endif
        END DO
      END DO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          snhflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rain)
# ifdef EMINUSP
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          evap)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          EminusP)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,isalt))
# endif
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rain, snhflx)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    LW_down, SW_down)
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zair)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pott, tair, qair)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    uwind, vwind)
#  ifdef SNOWFALL
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rain_lq, snow)
#  endif
#  ifdef WW3_MODEL
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Tauicex, Tauicey)
#  endif              
#  ifdef EMINUSP
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    evap, EminusP,                                &
     &                    stflx(:,:,isalt))
#  endif
# endif
      RETURN
      END SUBROUTINE atm2oi_flux_tile
#endif
      END MODULE atm2oi_flux_mod
