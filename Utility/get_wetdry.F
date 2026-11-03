#include "cppdefs.h"
#ifdef WET_DRY
      SUBROUTINE get_wetdry (ng, model, IniRec)
!
!svn $Id: get_wetdry.F 927 2018-10-16 03:51:56Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  If restart, reads wetting and drying masks from restart (initial)   !
!  NetCDF file.
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
#ifdef WEC_VF
      USE mod_ocean
#endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE exchange_2d_mod
#ifdef WEC_VF
      USE exchange_3d_mod
      USE mp_exchange_mod, ONLY : mp_exchange3d
      USE nf_fread3d_mod,  ONLY : nf_fread3d
#endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      USE nf_fread2d_mod,  ONLY : nf_fread2d
      USE strings_mod,     ONLY : find_string
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, IniRec
!
!  Local variable declarations.
!
      integer :: tile, LBi, UBi, LBj, UBj
      integer :: gtype, i, j,status
      integer :: varid_pmask, varid_rmask, varid_umask, varid_vmask,    &
     &           varid_umaskd, varid_vmaskd, varid_zetat,varid_ustokes, &
     &           varid_vstokes
      integer :: Vsize(4)

      real(dp), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin

      character (len=256) :: ncname
!
      SourceFile=__FILE__
!
!-----------------------------------------------------------------------
!  Inquire about the contents of restart NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
      ncname=INI(ng)%name
!
!  Open restart (initial) NetCDF for read/write.
!
      IF (INI(ng)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 1, INI(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Check grid file dimensions for consitency
!
      CALL netcdf_check_dim (ng, model, ncname,                         &
     &                       ncid = INI(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = INI(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idPwet)),        &
     &                     varid_pmask)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idPwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idRwet)),        &
     &                     varid_rmask)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idRwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idUwet)),        &
     &                     varid_umask)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idUwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idVwet)),        &
     &                     varid_vmask)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idVwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idUwetd)),       &
     &                     varid_umaskd)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idUwetd)),           &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF    
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idVwetd)),       &
     &                     varid_vmaskd)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idVwetd)),           &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#ifdef WEC_VF
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idWztt)),       &
     &                     varid_zetat)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idWztt)),           &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idU3Sd)),       &
     &                     varid_ustokes)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idU3Sd)),           &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,TRIM(Vname(1,idV3Sd)),       &
     &                     varid_vstokes)) THEN
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idV3Sd)),           &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF

#endif    
!-----------------------------------------------------------------------
!  Read in wet/dry mask from rstart file.
!-----------------------------------------------------------------------
!
!  Set 2D arrays bounds.
!
# ifdef DISTRIBUTE
      tile=MyRank
# else
      tile=-1
# endif
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Read in wet/dry mask at PSI-points.
!
      gtype=p2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idPwet), varid_pmask,                   &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % pmask,                               &
     &                  GRID(ng) % pmask_wet)
      IF (FoundError(status, nf90_noerr, __LINE__,                      &
     &               __FILE__)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idPwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idPwet)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in wet/dry mask at RHO-points.
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idRwet), varid_rmask,                   &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % rmask,                               &
     &                  GRID(ng) % rmask_wet)
      IF (FoundError(status, nf90_noerr, __LINE__,                      &
     &               __FILE__)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idRwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idRwet)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in wet/dry mask at U-points.
!
      gtype=u2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idUwet), varid_umask,                   &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % umask,                               &
     &                  GRID(ng) % umask_wet)
      IF (FoundError(status, nf90_noerr, __LINE__,                      &
     &               __FILE__)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idUwet)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in wet/dry mask at V-points.
!
      gtype=v2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idVwet), varid_vmask,                   &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % vmask,                               &
     &                  GRID(ng) % vmask_wet)
      IF (FoundError(status, nf90_noerr, __LINE__,                      &
     &               __FILE__)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVwet)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idVwet)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in wet/dry diffusion coefficient at U-points.
!
      gtype=u2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idUwetd), varid_umaskd,                 &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % umask,                               &
     &                  GRID(ng) % umask_diff)
      IF (FoundError(status, nf90_noerr, __LINE__,                      &
     &               __FILE__)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUwetd)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idUwetd)), ng, TRIM(ncname),   &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in wet/dry diffusion coefficient at V-points.
!
      gtype=v2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idVwetd), varid_vmaskd,                 &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % vmask,                               &
     &                  GRID(ng) % vmask_diff)
      IF (FoundError(status, nf90_noerr, __LINE__,                      &
     &               __FILE__)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVwetd)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idVwetd)), ng, TRIM(ncname),   &
     &                      Fmin, Fmax
        END IF
      END IF      
#ifdef WEC_VF
!
!  Read in zetat at Rho-points.
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idWztt), varid_zetat,                   &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  OCEAN(ng) % zetat)
      IF (FoundError(status, nf90_noerr, __LINE__,                      &
     &               __FILE__)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idWztt)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idWztt)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in u_stokes at U-points.
!
      gtype=u3dvar
      status=nf_fread3d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idU3Sd), varid_ustokes,                 &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % umask,                               &
# endif
     &                  OCEAN(ng) % u_stokes)
      IF (FoundError(status, nf90_noerr, __LINE__,                      &
     &               __FILE__)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idU3Sd)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idU3Sd)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
!
!  Read in v_stokes at V-points.
!
      gtype=v3dvar
      status=nf_fread3d(ng, model, ncname, INI(ng)%ncid,                &
     &                  Vname(1,idV3Sd), varid_vstokes,                 &
     &                  IniRec, gtype, Vsize,                           &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  OCEAN(ng) % v_stokes)
      IF (FoundError(status, nf90_noerr, __LINE__,                      &
     &               __FILE__)) THEN
        IF (Master) WRITE (stdout,30) TRIM(Vname(1,idV3Sd)),            &
     &                                TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,40) TRIM(Vname(2,idV3Sd)), ng, TRIM(ncname),    &
     &                      Fmin, Fmax
        END IF
      END IF
#endif          
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % pmask_wet)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % rmask_wet)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % umask_wet)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % vmask_wet)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % umask_diff)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % vmask_diff)
#ifdef WEC_VF
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(ng) % zetat)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          OCEAN(ng) % u_stokes)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          OCEAN(ng) % v_stokes)
#endif
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    GRID(ng) % pmask_wet,                         &
     &                    GRID(ng) % rmask_wet,                         &
     &                    GRID(ng) % umask_wet,                         &
     &                    GRID(ng) % vmask_wet)
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    GRID(ng) % umask_diff,                        &
     &                    GRID(ng) % vmask_diff)                        
#  ifdef WEC_VF
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    OCEAN(ng) % zetat)
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    OCEAN(ng) % u_stokes,                         &
     &                    OCEAN(ng) % v_stokes)                 
#  endif
# endif
!
  10  FORMAT (/,' GET_WETDRY - unable to open grid NetCDF file: ',a)
  20  FORMAT (/,' GET_WETDRY - unable to find grid variable: ',a,       &
     &        /,14x,'in grid NetCDF file: ',a)
  30  FORMAT (/,' GET_WETDRY - error while reading variable: ',a,       &
     &        /,12x,'in grid NetCDF file: ',a)
  40  FORMAT (3x,' GET_WETDRY  - ',a,/,19x,                             &
     &        '(Grid = ',i2.2,', File: ',a,')',/,19x,                   &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
#else
      SUBROUTINE get_wetdry
#endif
      RETURN
      END SUBROUTINE get_wetdry
