!
! Copyright (C) 2008-2012 Georgy Samsonidze
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Converts the output files produced by pw.x to the input files for BerkeleyGW.
!
!-------------------------------------------------------------------------------
!
! BerkeleyGW, Copyright (c) 2011, The Regents of the University of
! California, through Lawrence Berkeley National Laboratory (subject to
! receipt of any required approvals from the U.S. Dept. of Energy).
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
! (1) Redistributions of source code must retain the above copyright
! notice, this list of conditions and the following disclaimer.
!
! (2) Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! (3) Neither the name of the University of California, Lawrence
! Berkeley National Laboratory, U.S. Dept. of Energy nor the names of
! its contributors may be used to endorse or promote products derived
! from this software without specific prior written permission.
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
!
! You are under no obligation whatsoever to provide any bug fixes,
! patches, or upgrades to the features, functionality or performance of
! the source code ("Enhancements") to anyone; however, if you choose to
! make your Enhancements available either publicly, or directly to
! Lawrence Berkeley National Laboratory, without imposing a separate
! written license agreement for such Enhancements, then you hereby grant
! the following license: a  non-exclusive, royalty-free perpetual
! license to install, use, modify, prepare derivative works, incorporate
! into other computer software, distribute, and sublicense such
! enhancements or derivative works thereof, in binary and source code
! form.
!
!-------------------------------------------------------------------------------
!
! pw2bgw subroutines:
!
! write_wfng  - generates complex wavefunctions in G-space (normalized to 1)
! real_wfng   - constructs real wavefunctions by applying the Gram-Schmidt
!               process (called from write_wfng)
! write_rhog  - generates real/complex charge density in G-space
!               (units of the number of electronic states per unit cell)
! calc_rhog   - computes charge density by summing over a subset of occupied
!               bands (called from write_rhog), destroys charge density
! write_vxcg  - generates real/complex exchange-correlation potential in
!               G-space (units of Rydberg) [only local part of Vxc]
! write_vxc0  - prints real/complex exchange-correlation potential at G=0
!               (units of eV) [only local part of Vxc]
! write_vxc_r - calculates matrix elements of exchange-correlation potential
!               in R-space (units of eV) [only local part of Vxc]
! write_vxc_g - calculates matrix elements of exchange-correlation potential
!               in G-space (units of eV) [supports non-local Vxc]
! write_vscg  - generates real/complex self-consistent potential in G-space
!               (units of Rydberg) [only local part of Vsc]
! write_vkbg  - generates complex Kleinman-Bylander projectors in G-space
!               (units of Rydberg)
! check_inversion - checks whether real/complex version is appropriate
!               (called from everywhere)
!
! Quantum ESPRESSO stores the wavefunctions in is-ik-ib-ig order
! BerkeleyGW stores the wavefunctions in ik-ib-is-ig order
! the outer loop is over is(QE)/ik(BGW) and the inner loop is over ig
! ik = k-point index, is = spin index, ib = band index, ig = G-vector index
!
! write_wfng reverts the order of is and ik using smap and kmap arrays,
! distributes wavefunctions over processors by ig (either real case or
! spin-polarized case), calls real_wfng that applies the Gram-Schmidt
! process (real case), reverts the order of is and ib (spin-polarized
! case), and writes wavefunctions to disk
!
!-------------------------------------------------------------------------------

PROGRAM pw2bgw

  USE constants, ONLY : eps12
  USE control_flags, ONLY : gamma_only
  USE environment, ONLY : environment_start, environment_end
  USE io_files, ONLY : prefix, tmp_dir
  !USE io_global, ONLY : ionode, ionode_id        !FZ: commented
  USE io_global,        ONLY : ionode, ionode_id, stdout    !FZ: test if need this when output to pp_out
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_bcast
  USE mp_world, ONLY : world_comm
  USE mp_global, ONLY : mp_startup
  USE paw_variables, ONLY : okpaw
  USE scf, ONLY : rho_core, rhog_core
  USE uspp, ONLY : okvan

  IMPLICIT NONE

  character(len=6) :: codename = 'PW2BGW'

  integer :: real_or_complex
  character ( len = 9 ) :: symm_type
  logical :: wfng_flag
  character ( len = 256 ) :: wfng_file
  logical :: wfng_kgrid
  integer :: wfng_nk1
  integer :: wfng_nk2
  integer :: wfng_nk3
  real (DP) :: wfng_dk1
  real (DP) :: wfng_dk2
  real (DP) :: wfng_dk3
  logical :: wfng_occupation
  integer :: wfng_nvmin
  integer :: wfng_nvmax
  logical :: rhog_flag
  character ( len = 256 ) :: rhog_file
  integer :: rhog_nvmin
  integer :: rhog_nvmax
  logical :: vxcg_flag
  character ( len = 256 ) :: vxcg_file
  logical :: vxc0_flag
  character ( len = 256 ) :: vxc0_file
  logical :: vxc_flag
  character ( len = 256 ) :: vxc_file
  logical :: v_h_flag                 !FZ: for metaGGA output VHartree
  character ( len = 256 ) :: v_h_file !FZ: for metaGGA output VHartree
  logical :: ekin_flag                 !FZ: for metaGGA output kinetic energy
  character ( len = 256 ) :: ekin_file !FZ: for metaGGA output kinetic energy
  logical :: vltot_flag                 !FZ: for metaGGA output local potential energy
  character ( len = 256 ) :: vltot_file !FZ: for metaGGA output local potential energy
  logical :: vnl_flag                 !FZ: for metaGGA output nonlocal potential energy
  character ( len = 256 ) :: vnl_file !FZ: for metaGGA output nonlocal potential energy
  character :: vxc_integral
  integer :: vxc_diag_nmin
  integer :: vxc_diag_nmax
  integer :: vxc_offdiag_nmin
  integer :: vxc_offdiag_nmax
  integer :: v_h_diag_nmin        !FZ: for metaGGA output VHartree
  integer :: v_h_diag_nmax        !FZ: for metaGGA output VHartree
  integer :: v_h_offdiag_nmin     !FZ: for metaGGA output VHartree
  integer :: v_h_offdiag_nmax     !FZ: for metaGGA output VHartree
  integer :: ekin_diag_nmin        !FZ: for metaGGA output kinetic energy
  integer :: ekin_diag_nmax        !FZ: for metaGGA output kinetic energy
  integer :: vltot_diag_nmin        !FZ: for metaGGA output local potential energy
  integer :: vltot_diag_nmax        !FZ: for metaGGA output local potential energy
  integer :: vltot_offdiag_nmin     !FZ: for metaGGA output local potential energy
  integer :: vltot_offdiag_nmax     !FZ: for metaGGA output local potential energy
  integer :: vnl_diag_nmin        !FZ: for metaGGA output local potential energy
  integer :: vnl_diag_nmax        !FZ: for metaGGA output local potential energy
  integer :: vnl_offdiag_nmin     !FZ: for metaGGA output local potential energy
  integer :: vnl_offdiag_nmax     !FZ: for metaGGA output local potential energy
  logical :: vxc_zero_rho_core
  logical :: vscg_flag
  character ( len = 256 ) :: vscg_file
  logical :: vkbg_flag
  character ( len = 256 ) :: vkbg_file
  character ( len = 256 ) :: outdir

  NAMELIST / input_pw2bgw / prefix, outdir, &
    real_or_complex, symm_type, wfng_flag, wfng_file, wfng_kgrid, &
    wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3, &
    wfng_occupation, wfng_nvmin, wfng_nvmax, rhog_flag, rhog_file, &
    rhog_nvmin, rhog_nvmax, vxcg_flag, vxcg_file, vxc0_flag, vxc0_file, &
    vxc_flag, vxc_file, vxc_integral, vxc_diag_nmin, vxc_diag_nmax, &
    vxc_offdiag_nmin, vxc_offdiag_nmax, vxc_zero_rho_core, &
    vscg_flag, vscg_file, vkbg_flag, vkbg_file, v_h_flag, v_h_file, &   !FZ: for metaGGA output VHartree added v_h_flag, v_h_file
    ekin_flag, ekin_file, v_h_diag_nmin, v_h_diag_nmax, v_h_offdiag_nmin, v_h_offdiag_nmax, &  !FZ: for metaGGA output hartree
    ekin_diag_nmin, ekin_diag_nmax, vltot_flag, vltot_file, vltot_diag_nmin, vltot_diag_nmax, & !FZ: for metaGGA output local potential
    vltot_offdiag_nmin, vltot_offdiag_nmax, vnl_flag, vnl_file, vnl_diag_nmin, vnl_diag_nmax, &   !FZ: for metaGGA output local potential
    vnl_offdiag_nmin, vnl_offdiag_nmax  !FZ: for metaGGA output nonlocal potential

  integer :: ii, ios
  character ( len = 256 ) :: output_file_name

  character (len=256), external :: trimcheck
  character (len=1), external :: lowercase

#if defined(__MPI)
  CALL mp_startup ( )
#endif

  CALL environment_start ( codename )

  prefix = 'prefix'
  CALL get_environment_variable ( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM ( outdir ) == ' ' ) outdir = './'
  real_or_complex = 2
  symm_type = 'cubic'
  wfng_flag = .FALSE.
  wfng_file = 'WFN'
  wfng_kgrid = .FALSE.
  wfng_nk1 = 0
  wfng_nk2 = 0
  wfng_nk3 = 0
  wfng_dk1 = 0.0D0
  wfng_dk2 = 0.0D0
  wfng_dk3 = 0.0D0
  wfng_occupation = .FALSE.
  wfng_nvmin = 0
  wfng_nvmax = 0
  rhog_flag = .FALSE.
  rhog_file = 'RHO'
  rhog_nvmin = 0
  rhog_nvmax = 0
  vxcg_flag = .FALSE.
  vxcg_file = 'VXC'
  vxc0_flag = .FALSE.
  vxc0_file = 'vxc0.dat'
  vxc_flag = .FALSE.
  vxc_file = 'vxc.dat'
  v_h_flag = .TRUE.    !FZ: for metaGGA output VHartree 
  v_h_file = 'v_hartree.dat'  !FZ: for metaGGA output VHartree
  ekin_flag = .TRUE.    !FZ: for metaGGA output kinetic energy 
  ekin_file = 'ekin.dat'  !FZ: for metaGGA output kinetic energy
  vltot_flag = .TRUE.    !FZ: for metaGGA output local potential energy 
  vltot_file = 'vsc.dat'  !FZ: for metaGGA output local potential energy
  vnl_flag = .TRUE.    !FZ: for metaGGA output nonlocal potential energy 
  vnl_file = 'vnl.dat'  !FZ: for metaGGA output nonlocal potential energy
  vxc_integral = 'g'
  vxc_diag_nmin = 0
  vxc_diag_nmax = 0
  vxc_offdiag_nmin = 0
  vxc_offdiag_nmax = 0
  vxc_zero_rho_core = .TRUE.
  vscg_flag = .FALSE.
  vscg_file = 'VSC'
  vkbg_flag = .FALSE.
  vkbg_file = 'VKB'

  IF ( ionode ) THEN
    CALL input_from_file ( )
    READ ( 5, input_pw2bgw, iostat = ios )
    IF ( ios /= 0 ) CALL errore ( codename, 'input_pw2bgw', abs ( ios ) )

    DO ii = 1, LEN_TRIM (symm_type)
      symm_type(ii:ii) = lowercase (symm_type(ii:ii))
    END DO
    DO ii = 1, LEN_TRIM (vxc_integral)
      vxc_integral(ii:ii) = lowercase (vxc_integral(ii:ii))
    END DO

    IF ( real_or_complex /= 1 .AND. real_or_complex /= 2 ) &
      CALL errore ( codename, 'real_or_complex', 1 )
    IF ( symm_type /= 'cubic' .AND. symm_type /= 'hexagonal' ) &
      CALL errore ( codename, 'symm_type', 1 )
    IF ( vxc_integral /= 'r' .AND. vxc_integral /= 'g' ) &
      CALL errore ( codename, 'vxc_integral', 1 )
  ENDIF

  tmp_dir = trimcheck ( outdir )
  CALL mp_bcast ( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast ( prefix, ionode_id, world_comm )
  CALL mp_bcast ( real_or_complex, ionode_id, world_comm )
  CALL mp_bcast ( symm_type, ionode_id, world_comm )
  CALL mp_bcast ( wfng_flag, ionode_id, world_comm )
  CALL mp_bcast ( wfng_file, ionode_id, world_comm )
  CALL mp_bcast ( wfng_kgrid, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nk1, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nk2, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nk3, ionode_id, world_comm )
  CALL mp_bcast ( wfng_dk1, ionode_id, world_comm )
  CALL mp_bcast ( wfng_dk2, ionode_id, world_comm )
  CALL mp_bcast ( wfng_dk3, ionode_id, world_comm )
  CALL mp_bcast ( wfng_occupation, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nvmin, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nvmax, ionode_id, world_comm )
  CALL mp_bcast ( rhog_flag, ionode_id, world_comm )
  CALL mp_bcast ( rhog_file, ionode_id, world_comm )
  CALL mp_bcast ( rhog_nvmin, ionode_id, world_comm )
  CALL mp_bcast ( rhog_nvmax, ionode_id, world_comm )
  CALL mp_bcast ( vxcg_flag, ionode_id, world_comm )
  CALL mp_bcast ( vxcg_file, ionode_id, world_comm )
  CALL mp_bcast ( vxc0_flag, ionode_id, world_comm )
  CALL mp_bcast ( vxc0_file, ionode_id, world_comm )
  CALL mp_bcast ( vxc_flag, ionode_id, world_comm )
  CALL mp_bcast ( vxc_integral, ionode_id, world_comm )
  CALL mp_bcast ( vxc_file, ionode_id, world_comm )
  !CALL mp_bcast ( v_h_flag, ionode_id, world_comm )     !FZ: for metaGGA output VHartree
  !CALL mp_bcast ( v_h_file, ionode_id, world_comm )     !FZ: for metaGGA output VHartree
  !CALL mp_bcast ( ekin_flag, ionode_id, world_comm )     !FZ: for metaGGA output kinetic energy
  !CALL mp_bcast ( ekin_file, ionode_id, world_comm )     !FZ: for metaGGA output kinetic energy
  !CALL mp_bcast ( vltot_flag, ionode_id, world_comm )     !FZ: for metaGGA output local potential energy
  !CALL mp_bcast ( vltot_file, ionode_id, world_comm )     !FZ: for metaGGA output local potential energy
  !CALL mp_bcast ( vnl_flag, ionode_id, world_comm )     !FZ: for metaGGA output nonlocal potential energy
  !CALL mp_bcast ( vnl_file, ionode_id, world_comm )     !FZ: for metaGGA output nonlocal potential energy
  CALL mp_bcast ( vxc_diag_nmin, ionode_id, world_comm )
  CALL mp_bcast ( vxc_diag_nmax, ionode_id, world_comm )
  CALL mp_bcast ( vxc_offdiag_nmin, ionode_id, world_comm )
  CALL mp_bcast ( vxc_offdiag_nmax, ionode_id, world_comm )
  !CALL mp_bcast ( v_h_diag_nmin, ionode_id, world_comm )    !FZ: for metaGGA output VHartree
  !CALL mp_bcast ( v_h_diag_nmax, ionode_id, world_comm )    !FZ: for metaGGA output VHartree
  !CALL mp_bcast ( v_h_offdiag_nmin, ionode_id, world_comm )    !FZ: for metaGGA output VHartree
  !CALL mp_bcast ( v_h_offdiag_nmax, ionode_id, world_comm )    !FZ: for metaGGA output VHartree
  !CALL mp_bcast ( ekin_diag_nmin, ionode_id, world_comm )    !FZ: for metaGGA output kinetic energy
  !CALL mp_bcast ( ekin_diag_nmax, ionode_id, world_comm )    !FZ: for metaGGA output kinetic energy
  !CALL mp_bcast ( vltot_diag_nmin, ionode_id, world_comm )    !FZ: for metaGGA output local potential energy
  !CALL mp_bcast ( vltot_diag_nmax, ionode_id, world_comm )    !FZ: for metaGGA output local potential energy
  !CALL mp_bcast ( vltot_offdiag_nmin, ionode_id, world_comm )  !FZ: for metaGGA output local potential energy
  !CALL mp_bcast ( vltot_offdiag_nmax, ionode_id, world_comm )  !FZ: for metaGGA output local potential energy
  !CALL mp_bcast ( vnl_diag_nmin, ionode_id, world_comm )    !FZ: for metaGGA output nonlocal potential energy
  !CALL mp_bcast ( vnl_diag_nmax, ionode_id, world_comm )    !FZ: for metaGGA output nonlocal potential energy
  !CALL mp_bcast ( vnl_offdiag_nmin, ionode_id, world_comm )  !FZ: for metaGGA output nonlocal potential energy
  !CALL mp_bcast ( vnl_offdiag_nmax, ionode_id, world_comm )  !FZ: for metaGGA output nonlocal potential energy
  CALL mp_bcast ( vxc_zero_rho_core, ionode_id, world_comm )
  CALL mp_bcast ( vscg_flag, ionode_id, world_comm )
  CALL mp_bcast ( vscg_file, ionode_id, world_comm )
  CALL mp_bcast ( vkbg_flag, ionode_id, world_comm )
  CALL mp_bcast ( vkbg_file, ionode_id, world_comm )

  CALL read_file ( )

  v_h_diag_nmin = vxc_diag_nmin       !FZ: for metaGGA output VHartree
  v_h_diag_nmax = vxc_diag_nmax       !FZ: for metaGGA output VHartree
  v_h_offdiag_nmin = vxc_offdiag_nmin    !FZ: for metaGGA output VHartree
  v_h_offdiag_nmax = vxc_offdiag_nmax    !FZ: for metaGGA output VHartree
  ekin_diag_nmin = vxc_diag_nmin       !FZ: for metaGGA output kinetic energy
  ekin_diag_nmax = vxc_diag_nmax       !FZ: for metaGGA output kinetic energy
  vltot_diag_nmin = vxc_diag_nmin       !FZ: for metaGGA output local potential energy
  vltot_diag_nmax = vxc_diag_nmax       !FZ: for metaGGA output local potential energy
  vltot_offdiag_nmin = vxc_offdiag_nmin    !FZ: for metaGGA output local potential energy
  vltot_offdiag_nmax = vxc_offdiag_nmax    !FZ: for metaGGA output local potential energy
  vnl_diag_nmin = vxc_diag_nmin       !FZ: for metaGGA output local potential energy
  vnl_diag_nmax = vxc_diag_nmax       !FZ: for metaGGA output local potential energy
  vnl_offdiag_nmin = vxc_offdiag_nmin    !FZ: for metaGGA output local potential energy
  vnl_offdiag_nmax = vxc_offdiag_nmax    !FZ: for metaGGA output local potential energy

  if (ionode) then
    if (MAX (MAXVAL (ABS (rho_core (:) ) ), MAXVAL (ABS (rhog_core (:) ) ) ) &
      .LT. eps12) then
      WRITE ( 6, '(/,5x,"NLCC is absent")' )
    else
      WRITE ( 6, '(/,5x,"NLCC is present")' )
    endif
  endif
  if (okvan) call errore ( 'pw2bgw', 'BGW cannot use USPP.', 3 )
  if (okpaw) call errore ( 'pw2bgw', 'BGW cannot use PAW.', 4 )
  if (gamma_only) call errore ( 'pw2bgw', 'BGW cannot use gamma-only run.', 5 )
  if (nspin == 4) call errore ( 'pw2bgw', 'BGW cannot use spinors.', 6 )
  if (real_or_complex == 1 .AND. vxc_flag .AND. vxc_offdiag_nmax > 0) &
    call errore ( 'pw2bgw', 'Off-diagonal matrix elements of Vxc ' // &
    'with real wavefunctions are not implemented, compute them in ' // &
    'Sigma using VXC.', 7)

  CALL openfil_pp ( )

  if ( ionode ) WRITE ( 6, '("")' )

  IF ( ionode ) WRITE (stdout, '(" Test that stdout output something")')    !FZ: test
  IF ( ionode ) WRITE (6, '(" Test that WRITE(6, ...) output something")')  !FZ: test

  IF ( wfng_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( wfng_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_wfng")' )
    CALL start_clock ( 'write_wfng' )
    CALL write_wfng ( output_file_name, real_or_complex, symm_type, &
      wfng_kgrid, wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, &
      wfng_dk3, wfng_occupation, wfng_nvmin, wfng_nvmax )
    CALL stop_clock ( 'write_wfng' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_wfng",/)' )
  ENDIF

  IF ( vxcg_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( vxcg_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vxcg")' )
    CALL start_clock ( 'write_vxcg' )
    CALL write_vxcg ( output_file_name, real_or_complex, symm_type, &
      vxc_zero_rho_core )
    CALL stop_clock ( 'write_vxcg' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vxcg",/)' )
  ENDIF

  IF ( vxc0_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( vxc0_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc0")' )
    CALL start_clock ( 'write_vxc0' )
    CALL write_vxc0 ( output_file_name, vxc_zero_rho_core )
    CALL stop_clock ( 'write_vxc0' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc0",/)' )
  ENDIF

  IF ( vxc_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( vxc_file )
    IF ( ionode ) WRITE (6, '(" Test that write_vxc_r called")')  !FZ: test
    IF ( vxc_integral .EQ. 'r' ) THEN
      IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc_r")' )
      CALL start_clock ( 'write_vxc_r' )
      CALL write_vxc_r ( output_file_name, &
        vxc_diag_nmin, vxc_diag_nmax, &
        vxc_offdiag_nmin, vxc_offdiag_nmax, &
        vxc_zero_rho_core )
      CALL stop_clock ( 'write_vxc_r' )
      IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc_r",/)' )
    ENDIF
    IF ( vxc_integral .EQ. 'g' ) THEN
      !IF ( ionode ) WRITE (6, *) "vxc_diag_nmin = ", vxc_diag_nmin, "vxc_diag_nmax = ", vxc_diag_nmax, "vxc_offdiag_nmin = ", vxc_offdiag_nmin, "vxc_offdiag_nmax = ", vxc_offdiag_nmax, "vxc_zero_rho_core = ", vxc_zero_rho_core !FZ: test
      IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc_g")' )
      CALL start_clock ( 'write_vxc_g' )
      CALL write_vxc_g ( output_file_name, &
        vxc_diag_nmin, vxc_diag_nmax, &
        vxc_offdiag_nmin, vxc_offdiag_nmax, &
        vxc_zero_rho_core )
      CALL stop_clock ( 'write_vxc_g' )
      IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc_g",/)' )
    ENDIF
  ENDIF
  
  !FZ:   This block is for metaGGA to output Hartree energy for each bands
  IF ( v_h_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( v_h_file )
    IF ( vxc_integral .EQ. 'r' ) THEN              !FZ this subblock is reserved for v_h_r function if needed
      CALL errore ("pw2bgw", "write_v_h_r has vxc_integral .EQ. 'r', write_v_h_r has not been written", 1)
      !IF ( ionode ) WRITE ( 6, '(5x,"call write_v_h_r")' )
      !CALL start_clock ( 'write_v_h_r' )
      !CALL write_v_h_r ( output_file_name, &
      !  v_h_diag_nmin, v_h_diag_nmax, &
      !  v_h_offdiag_nmin, v_h_offdiag_nmax)
      !CALL stop_clock ( 'write_v_h_r' )
      !IF ( ionode ) WRITE ( 6, '(5x,"done write_v_h_r",/)' )
    ENDIF
    IF ( vxc_integral .EQ. 'g' ) THEN
      IF ( ionode ) WRITE (6, *) "v_h_diag_nmin = ", v_h_diag_nmin, "v_h_diag_nmax = ", v_h_diag_nmax, "v_h_offdiag_nmin = ", v_h_offdiag_nmin, "v_h_offdiag_nmax = ", v_h_offdiag_nmax !FZ: test
      IF ( ionode ) WRITE ( 6, '(5x,"call write_v_h_g")' )
      CALL start_clock ( 'write_v_h_g' )
      CALL write_v_h_g ( output_file_name, &
        v_h_diag_nmin, v_h_diag_nmax, &
        v_h_offdiag_nmin, v_h_offdiag_nmax)
      CALL stop_clock ( 'write_v_h_g' )
      IF ( ionode ) WRITE ( 6, '(5x,"done write_v_h_g",/)' )
    ENDIF
  ENDIF
  
  !FZ:   This block is for metaGGA to output local potential energy for each bands
  IF ( vltot_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( vltot_file )
    IF ( vxc_integral .EQ. 'r' ) THEN              !FZ this subblock is reserved for v_h_r function if needed
      CALL errore ("pw2bgw", "write_vltot has vxc_integral .EQ. 'r', write_vltot_r has not been written", 1)
      !IF ( ionode ) WRITE ( 6, '(5x,"call write_v_h_r")' )
      !CALL start_clock ( 'write_v_h_r' )
      !CALL write_v_h_r ( output_file_name, &
      !  v_h_diag_nmin, v_h_diag_nmax, &
      !  v_h_offdiag_nmin, v_h_offdiag_nmax)
      !CALL stop_clock ( 'write_v_h_r' )
      !IF ( ionode ) WRITE ( 6, '(5x,"done write_v_h_r",/)' )
    ENDIF
    IF ( vxc_integral .EQ. 'g' ) THEN
      IF ( ionode ) WRITE (6, *) "vltot_diag_nmin = ", vltot_diag_nmin, "vltot_diag_nmax = ", vltot_diag_nmax, "vltot_offdiag_nmin = ", vltot_offdiag_nmin, "vltot_offdiag_nmax = ", vltot_offdiag_nmax !FZ: test
      IF ( ionode ) WRITE ( 6, '(5x,"call write_vltot")' )
      CALL start_clock ( 'write_vltot' )
      CALL write_vltot ( output_file_name, &
        vltot_diag_nmin, vltot_diag_nmax, &
        vltot_offdiag_nmin, vltot_offdiag_nmax)
      CALL stop_clock ( 'write_vltot' )
      IF ( ionode ) WRITE ( 6, '(5x,"done write_vltot",/)' )
    ENDIF
  ENDIF
  
  !FZ:   This block is for metaGGA to output kinetic energy for each bands
  IF ( ekin_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( ekin_file )
    !IF ( vxc_integral .EQ. 'g' ) THEN
      IF ( ionode ) WRITE (6, *) "ekin_diag_nmin = ", ekin_diag_nmin, "ekin_diag_nmax = ", ekin_diag_nmax !FZ: test
      IF ( ionode ) WRITE ( 6, '(5x,"call write_ekin_g")' )
      CALL start_clock ( 'write_ekin_g' )
      CALL write_ekin_g ( output_file_name, &
        ekin_diag_nmin, ekin_diag_nmax)
      CALL stop_clock ( 'write_ekin_g' )
      IF ( ionode ) WRITE ( 6, '(5x,"done write_ekin_g",/)' )
    !ENDIF
  ENDIF

  !FZ:   This block is for metaGGA to output nonlocal potential energy for each bands
  IF ( vnl_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( vnl_file )
    IF ( vxc_integral .EQ. 'r' ) THEN              !FZ this subblock is reserved for v_h_r function if needed
      CALL errore ("pw2bgw", "write_vnl has vxc_integral .EQ. 'r', write_vnl_r has not been written", 1)
      !IF ( ionode ) WRITE ( 6, '(5x,"call write_v_h_r")' )
      !CALL start_clock ( 'write_v_h_r' )
      !CALL write_v_h_r ( output_file_name, &
      !  v_h_diag_nmin, v_h_diag_nmax, &
      !  v_h_offdiag_nmin, v_h_offdiag_nmax)
      !CALL stop_clock ( 'write_v_h_r' )
      !IF ( ionode ) WRITE ( 6, '(5x,"done write_v_h_r",/)' )
    ENDIF
    IF ( vxc_integral .EQ. 'g' ) THEN
      IF ( ionode ) WRITE ( 6, '(5x,"call write_vnl")' )
      CALL start_clock ( 'write_vnl' )
      CALL write_vnl ( output_file_name, &
        vnl_diag_nmin, vnl_diag_nmax, &
        vnl_offdiag_nmin, vnl_offdiag_nmax)
      CALL stop_clock ( 'write_vnl' )
      IF ( ionode ) WRITE ( 6, '(5x,"done write_vnl",/)' )
    ENDIF
  ENDIF

  IF ( vscg_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( vscg_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vscg")' )
    CALL start_clock ( 'write_vscg' )
    CALL write_vscg ( output_file_name, real_or_complex, symm_type )
    CALL stop_clock ( 'write_vscg' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vscg",/)' )
  ENDIF

  IF ( vkbg_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( vkbg_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vkbg")' )
    CALL start_clock ( 'write_vkbg' )
    CALL write_vkbg ( output_file_name, symm_type, wfng_kgrid, wfng_nk1, &
      wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3 )
    CALL stop_clock ( 'write_vkbg' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vkbg",/)' )
  ENDIF

  ! since calc_rhog (called from write_rhog) destroys charge density,
  ! it must be called after v_xc (called from write_vxcg, write_vxc0,
  ! write_vxc_r, write_vxc_g)
  IF ( rhog_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( rhog_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_rhog")' )
    CALL start_clock ( 'write_rhog' )
    CALL write_rhog ( output_file_name, real_or_complex, symm_type, &
      rhog_nvmin, rhog_nvmax )
    CALL stop_clock ( 'write_rhog' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_rhog",/)' )
  ENDIF

  IF ( ionode ) WRITE ( 6, * )
  IF ( wfng_flag ) CALL print_clock ( 'write_wfng' )
  IF ( rhog_flag ) CALL print_clock ( 'write_rhog' )
  IF ( vxcg_flag ) CALL print_clock ( 'write_vxcg' )
  IF ( vxc0_flag ) CALL print_clock ( 'write_vxc0' )
  IF ( vxc_flag ) THEN
    IF ( vxc_integral .EQ. 'r' ) CALL print_clock ( 'write_vxc_r' )
    IF ( vxc_integral .EQ. 'g' ) CALL print_clock ( 'write_vxc_g' )
  ENDIF
  IF ( vscg_flag ) CALL print_clock ( 'write_vscg' )
  IF ( vkbg_flag ) CALL print_clock ( 'write_vkbg' )
  IF ( wfng_flag .AND. real_or_complex .EQ. 1 ) THEN
    IF ( ionode ) WRITE ( 6, '(/,5x,"Called by write_wfng:")' )
    CALL print_clock ( 'real_wfng' )
  ENDIF

  CALL environment_end ( codename )

  CALL stop_pp ( )

  ! this is needed because openfil is called above
  CALL close_files ( .false. )

  STOP

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE write_wfng ( output_file_name, real_or_complex, symm_type, &
  wfng_kgrid, wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, &
  wfng_dk3, wfng_occupation, wfng_nvmin, wfng_nvmax )

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  USE constants, ONLY : pi, tpi, eps6
  USE fft_base, ONLY : dfftp
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, g, mill, ecutrho
  USE io_files, ONLY : iunwfc, nwordwfc
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat, atm, ityp, tau
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, wk, ngk, nks, nkstot, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum, mp_max, mp_get, mp_bcast, mp_barrier
  USE mp_pools, ONLY : me_pool, root_pool, npool, nproc_pool, &
    intra_pool_comm, inter_pool_comm
  USE mp_wave, ONLY : mergewf
  USE mp_world, ONLY : mpime, nproc, world_comm
  USE mp_bands, ONLY : intra_bgrp_comm, nbgrp
  USE start_k, ONLY : nk1, nk2, nk3, k1, k2, k3
  USE symm_base, ONLY : s, ftau, nsym
  USE wavefunctions, ONLY : evc
  USE wvfct, ONLY : npwx, nbnd, et, wg
  USE gvecw, ONLY : ecutwfc
  USE matrix_inversion
#if defined(__MPI)
  USE parallel_include, ONLY : MPI_DOUBLE_COMPLEX
#endif

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  integer, intent (in) :: real_or_complex
  character ( len = 9 ), intent (in) :: symm_type
  logical, intent (in) :: wfng_kgrid
  integer, intent (in) :: wfng_nk1
  integer, intent (in) :: wfng_nk2
  integer, intent (in) :: wfng_nk3
  real (DP), intent (in) :: wfng_dk1
  real (DP), intent (in) :: wfng_dk2
  real (DP), intent (in) :: wfng_dk3
  logical, intent (in) :: wfng_occupation
  integer, intent (in) :: wfng_nvmin
  integer, intent (in) :: wfng_nvmax

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  logical :: proc_wf, bad_kgrid
  integer :: unit, i, j, k, cell_symmetry, nrecord
  integer :: id, ib, ik, iks, ike, is, ig, ierr
  integer :: nd, ntran, nb, nk_l, nk_g, ns, ng_l, ng_g
  integer :: npw, ngg, npw_g, npwx_g
  integer :: local_pw, ipsour, igwx, ngkdist_g, ngkdist_l
  real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: kmap ( : )
  integer, allocatable :: smap ( : )
  integer, allocatable :: ifmin ( : )
  integer, allocatable :: ifmax ( : )
  integer, allocatable :: itmp ( : )
  integer, allocatable :: ngk_g ( : )
  integer, allocatable :: ipmask ( : )
  integer, allocatable :: igwk ( : )
  integer, allocatable :: igwf_l2g ( : )
  integer, allocatable :: g_g ( :, : )
  integer, allocatable :: igk_l2g ( :, : )
  real (DP), allocatable :: et_g ( :, : )
  real (DP), allocatable :: wg_g ( :, : )
  real (DP), allocatable :: energy ( :, : )
  complex (DP), allocatable :: wfng ( : )
  complex (DP), allocatable :: wfng_buf ( :, : )
  complex (DP), allocatable :: wfng_dist ( :, :, : )

  INTEGER, EXTERNAL :: atomic_number, global_kpoint_index

  IF ( real_or_complex .EQ. 1 .OR. nspin .GT. 1 ) THEN
    proc_wf = .TRUE.
  ELSE
    proc_wf = .FALSE.
  ENDIF

  bad_kgrid = .FALSE.
  IF ( wfng_kgrid ) THEN
    IF ( wfng_nk1 .LE. 0 .OR. wfng_nk2 .LE. 0 .OR. wfng_nk3 .LE. 0 ) &
      bad_kgrid = .TRUE.
  ELSE
    IF ( nk1 .LE. 0 .OR. nk2 .LE. 0 .OR. nk3 .LE. 0 ) &
      bad_kgrid = .TRUE.
  ENDIF
  IF ( bad_kgrid .AND. ionode ) THEN
    WRITE ( 6, 101 )
  ENDIF

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  IF ( real_or_complex .EQ. 1 ) THEN
    WRITE ( stitle, '("WFN-Real",24X)' )
  ELSE
    WRITE ( stitle, '("WFN-Complex",21X)' )
  ENDIF

  unit = 4
  nrecord = 1
  nd = 3

  nb = nbnd
  nk_l = nks
  nk_g = nkstot
  ns = nspin
  ng_l = ngm
  ng_g = ngm_g

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  ALLOCATE ( kmap ( nk_g ) )
  ALLOCATE ( smap ( nk_g ) )

  DO i = 1, nk_g
    j = ( i - 1 ) / ns
    k = i - 1 - j * ns
    kmap ( i ) = j + k * ( nk_g / ns ) + 1
    smap ( i ) = k + 1
  ENDDO
  ierr = 0
  DO i = 1, nk_g
    ik = kmap ( i )
    is = smap ( i )
    IF ( ik .GE. iks .AND. ik .LE. ike .AND. is .NE. isk ( ik ) ) &
      ierr = ierr + 1
  ENDDO
  CALL mp_max ( ierr, world_comm )
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_wfng', 'smap', ierr )

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_wfng', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2 )
    t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( dfftp%nr1 )
    t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( dfftp%nr2 )
    t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( dfftp%nr3 )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

  ALLOCATE ( et_g ( nb, nk_g ) )

  DO ik = 1, nk_l
    DO ib = 1, nb
      et_g ( ib, ik ) = et ( ib, ik )
    ENDDO
  ENDDO
#if defined(__MPI)
  CALL poolrecover ( et_g, nb, nk_g, nk_l )
  CALL mp_bcast ( et_g, ionode_id, world_comm )
#endif

  ALLOCATE ( wg_g ( nb, nk_g ) )
  ALLOCATE ( ifmin ( nk_g ) )
  ALLOCATE ( ifmax ( nk_g ) )

  IF ( wfng_occupation ) THEN

    DO ik = 1, nk_g
      DO ib = 1, nb
        IF ( ib .GE. wfng_nvmin .AND. ib .LE. wfng_nvmax ) THEN
          wg_g ( ib, ik ) = 1.0D0
        ELSE
          wg_g ( ib, ik ) = 0.0D0
        ENDIF
      ENDDO
    ENDDO
    DO ik = 1, nk_g
      ifmin ( ik ) = wfng_nvmin
    ENDDO
    DO ik = 1, nk_g
      ifmax ( ik ) = wfng_nvmax
    ENDDO

  ELSE

    DO ik = 1, nk_l
      DO ib = 1, nb
        IF ( wk(ik) == 0.D0 ) THEN
          wg_g(ib,ik) = wg(ib,ik)
        ELSE
          wg_g(ib,ik) = wg(ib,ik) / wk(ik)
        ENDIF
      ENDDO
    ENDDO
#if defined(__MPI)
    CALL poolrecover ( wg_g, nb, nk_g, nk_l )
#endif
    DO ik = 1, nk_g
      ifmin ( ik ) = 0
    ENDDO
    DO ik = 1, nk_g
      ifmax ( ik ) = 0
    ENDDO
    DO ik = 1, nk_g
      DO ib = 1, nb
        IF ( wg_g( ib, ik ) .GT. 0.5D0 ) THEN
          IF ( ifmin ( ik ) .EQ. 0 ) ifmin ( ik ) = ib
          ifmax ( ik ) = ib
        ENDIF
      ENDDO
    ENDDO

  ENDIF

  ALLOCATE ( g_g ( nd, ng_g ) )

  DO ig = 1, ng_g
    DO id = 1, nd
      g_g ( id, ig ) = 0
    ENDDO
  ENDDO
  DO ig = 1, ng_l
    g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  CALL mp_sum ( g_g, intra_bgrp_comm )

  ALLOCATE ( igk_l2g ( npwx, nk_l ) )

  DO ik = 1, nk_l
    npw = ngk ( ik )
    DO ig = 1, npw
      igk_l2g ( ig, ik ) = ig_l2g ( igk_k (ig, ik) )
    ENDDO
    DO ig = npw + 1, npwx
      igk_l2g ( ig, ik ) = 0
    ENDDO
  ENDDO

  ALLOCATE ( ngk_g ( nk_g ) )

  DO ik = 1, nk_g
    ngk_g ( ik ) = 0
  ENDDO
  DO ik = 1, nk_l
    ngk_g ( ik + iks - 1 ) = ngk ( ik )
  ENDDO
  CALL mp_sum( ngk_g, inter_pool_comm )
  CALL mp_sum( ngk_g, intra_pool_comm )
  ngk_g = ngk_g / nbgrp

  npw_g = MAXVAL ( igk_l2g ( :, : ) )
  CALL mp_max( npw_g, intra_pool_comm )
  CALL mp_max( npw_g, inter_pool_comm )

  npwx_g = MAXVAL ( ngk_g ( : ) )

  CALL cryst_to_cart ( nk_g / ns, xk, at, - 1 )

  IF ( ionode ) THEN
    OPEN ( unit = unit, file = TRIM ( output_file_name ), &
      form = 'unformatted', status = 'replace' )
    WRITE ( unit ) stitle, sdate, stime
    WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho, &
      nk_g / ns, nb, npwx_g, ecutwfc
    IF ( wfng_kgrid ) THEN
      WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, wfng_nk1, wfng_nk2, wfng_nk3, &
        wfng_dk1, wfng_dk2, wfng_dk3
    ELSE
      WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, nk1, nk2, nk3, &
        0.5D0 * dble ( k1 ), 0.5D0 * dble ( k2 ), 0.5D0 * dble ( k3 )
    ENDIF
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) ( ngk_g ( ik ), ik = 1, nk_g / ns )
    WRITE ( unit ) ( wk ( ik ) * dble ( ns ) / 2.0D0, ik = 1, nk_g / ns )
    WRITE ( unit ) ( ( xk ( id, ik ), id = 1, nd ), ik = 1, nk_g / ns )
    WRITE ( unit ) ( ifmin ( ik ), ik = 1, nk_g )
    WRITE ( unit ) ( ifmax ( ik ), ik = 1, nk_g )
    WRITE ( unit ) ( ( et_g ( ib, ik ), ib = 1, nb ), ik = 1, nk_g )
    WRITE ( unit ) ( ( wg_g ( ib, ik ), ib = 1, nb ), ik = 1, nk_g )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
  ENDIF

  CALL cryst_to_cart ( nk_g / ns, xk, bg, 1 )

  DEALLOCATE ( wg_g )
  DEALLOCATE ( ifmax )
  DEALLOCATE ( ifmin )

  ALLOCATE ( igwk ( npwx_g ) )

  IF ( proc_wf ) THEN
    IF ( MOD ( npwx_g, nproc ) .EQ. 0 ) THEN
      ngkdist_l = npwx_g / nproc
    ELSE
      ngkdist_l = npwx_g / nproc + 1
    ENDIF
    ngkdist_g = ngkdist_l * nproc
    IF ( real_or_complex .EQ. 1 ) &
    ALLOCATE ( energy ( nb, ns ) )
    ALLOCATE ( wfng_buf ( ngkdist_g, ns ) )
    ALLOCATE ( wfng_dist ( ngkdist_l, nb, ns ) )
  ENDIF

  DO i = 1, nk_g

    ik = kmap ( i )
    is = smap ( i )

    IF ( real_or_complex .EQ. 1 ) THEN
      DO ib = 1, nb
        energy ( ib, is ) = et_g ( ib, i )
      ENDDO
    ENDIF

    DO j = 1, npwx_g
      igwk ( j ) = 0
    ENDDO
    ALLOCATE ( itmp ( npw_g ) )
    DO j = 1, npw_g
      itmp ( j ) = 0
    ENDDO
    IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
      DO ig = 1, ngk ( ik - iks + 1 )
        itmp ( igk_l2g ( ig, ik - iks + 1 ) ) = igk_l2g ( ig, ik - iks + 1 )
      ENDDO
    ENDIF
    CALL mp_sum( itmp, intra_bgrp_comm )
    !XXX
    CALL mp_sum( itmp, inter_pool_comm )
    !XXX
    ngg = 0
    DO ig = 1, npw_g
      IF ( itmp ( ig ) .EQ. ig ) THEN
        ngg = ngg + 1
        igwk ( ngg ) = ig
      ENDIF
    ENDDO
    DEALLOCATE ( itmp )

    IF ( ionode ) THEN
      IF ( is .EQ. 1 ) THEN
        WRITE ( unit ) nrecord
        WRITE ( unit ) ngk_g ( ik )
        WRITE ( unit ) ( ( g_g ( id, igwk ( ig ) ), id = 1, nd ), &
          ig = 1, ngk_g ( ik ) )
      ENDIF
    ENDIF

    local_pw = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
      CALL davcio ( evc, 2*nwordwfc, iunwfc, ik - iks + 1, - 1 )
      local_pw = ngk ( ik - iks + 1 )
    ENDIF

    ALLOCATE ( igwf_l2g ( local_pw ) )

    DO ig = 1, local_pw
      igwf_l2g ( ig ) = 0
    ENDDO
    DO ig = 1, local_pw
      ngg = igk_l2g ( ig, ik - iks + 1 )
      DO j = 1, ngk_g ( ik )
        IF ( ngg .EQ. igwk ( j ) ) THEN
          igwf_l2g ( ig ) = j
          EXIT
        ENDIF
      ENDDO
    ENDDO

    ALLOCATE ( ipmask ( nproc ) )
    DO j = 1, nproc
      ipmask ( j ) = 0
    ENDDO
    ipsour = ionode_id
    IF ( npool .GT. 1 ) THEN
      IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
        IF ( me_pool .EQ. root_pool ) ipmask ( mpime + 1 ) = 1
      ENDIF
      CALL mp_sum ( ipmask, world_comm )
      DO j = 1, nproc
        IF ( ipmask ( j ) .EQ. 1 ) ipsour = j - 1
      ENDDO
    ENDIF
    DEALLOCATE ( ipmask )

    igwx = 0
    ierr = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) &
      igwx = MAXVAL ( igwf_l2g ( 1 : local_pw ) )
    CALL mp_max ( igwx, intra_pool_comm )
    IF ( ipsour .NE. ionode_id ) &
      CALL mp_get ( igwx, igwx, mpime, ionode_id, ipsour, 1, world_comm )
    ierr = 0
    IF ( ik .GE. iks .AND. ik .LE. ike .AND. igwx .NE. ngk_g ( ik ) ) &
      ierr = 1
    CALL mp_max ( ierr, world_comm )
    IF ( ierr .GT. 0 ) &
      CALL errore ( 'write_wfng', 'igwx ngk_g', ierr )

    ALLOCATE ( wfng ( MAX ( 1, igwx ) ) )

    DO ib = 1, nb

      DO j = 1, igwx
        wfng ( j ) = ( 0.0D0, 0.0D0 )
      ENDDO
      IF ( npool .GT. 1 ) THEN
        IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
          CALL mergewf ( evc ( :, ib ), wfng, local_pw, igwf_l2g, &
            me_pool, nproc_pool, root_pool, intra_pool_comm )
        ENDIF
        IF ( ipsour .NE. ionode_id ) THEN
          CALL mp_get ( wfng, wfng, mpime, ionode_id, ipsour, ib, &
            world_comm )
        ENDIF
      ELSE
        CALL mergewf ( evc ( :, ib ), wfng, local_pw, igwf_l2g, &
          mpime, nproc, ionode_id, world_comm )
      ENDIF

      IF ( proc_wf ) THEN
        DO ig = 1, igwx
          wfng_buf ( ig, is ) = wfng ( ig )
        ENDDO
        DO ig = igwx + 1, ngkdist_g
          wfng_buf ( ig, is ) = ( 0.0D0, 0.0D0 )
        ENDDO
#if defined(__MPI)
        CALL mp_barrier ( world_comm )
        CALL MPI_Scatter ( wfng_buf ( :, is ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
        wfng_dist ( :, ib, is ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
        ionode_id, world_comm, ierr )
        IF ( ierr .GT. 0 ) &
          CALL errore ( 'write_wfng', 'mpi_scatter', ierr )
#else
        DO ig = 1, ngkdist_g
          wfng_dist ( ig, ib, is ) = wfng_buf ( ig, is )
        ENDDO
#endif
      ELSE
        IF ( ionode ) THEN
          WRITE ( unit ) nrecord
          WRITE ( unit ) ngk_g ( ik )
          WRITE ( unit ) ( wfng ( ig ), ig = 1, igwx )
        ENDIF
      ENDIF

    ENDDO

    DEALLOCATE ( wfng )
    DEALLOCATE ( igwf_l2g )

    IF ( proc_wf .AND. is .EQ. ns ) THEN
      IF ( real_or_complex .EQ. 1 ) THEN
        CALL start_clock ( 'real_wfng' )
        CALL real_wfng ( ik, ngkdist_l, nb, ns, energy, wfng_dist )
        CALL stop_clock ( 'real_wfng' )
      ENDIF
      DO ib = 1, nb
        DO is = 1, ns
#if defined(__MPI)
          CALL mp_barrier ( world_comm )
          CALL MPI_Gather ( wfng_dist ( :, ib, is ), ngkdist_l, &
          MPI_DOUBLE_COMPLEX, wfng_buf ( :, is ), ngkdist_l, &
          MPI_DOUBLE_COMPLEX, ionode_id, world_comm, ierr )
          IF ( ierr .GT. 0 ) &
            CALL errore ( 'write_wfng', 'mpi_gather', ierr )
#else
          DO ig = 1, ngkdist_g
            wfng_buf ( ig, is ) = wfng_dist ( ig, ib, is )
          ENDDO
#endif
        ENDDO
        IF ( ionode ) THEN
          WRITE ( unit ) nrecord
          WRITE ( unit ) ngk_g ( ik )
          IF ( real_or_complex .EQ. 1 ) THEN
            WRITE ( unit ) ( ( dble ( wfng_buf ( ig, is ) ), &
              ig = 1, igwx ), is = 1, ns )
          ELSE
            WRITE ( unit ) ( ( wfng_buf ( ig, is ), &
              ig = 1, igwx ), is = 1, ns )
          ENDIF
        ENDIF
      ENDDO
    ENDIF

  ENDDO

  DEALLOCATE ( igwk )
  DEALLOCATE ( ngk_g )
  DEALLOCATE ( igk_l2g )
  DEALLOCATE ( et_g )

  IF ( proc_wf ) THEN
    IF ( real_or_complex .EQ. 1 ) &
    DEALLOCATE ( energy )
    DEALLOCATE ( wfng_buf )
    DEALLOCATE ( wfng_dist )
  ENDIF

  IF ( ionode ) THEN
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( g_g )
  DEALLOCATE ( smap )
  DEALLOCATE ( kmap )

  CALL mp_barrier ( world_comm )

  RETURN

101 FORMAT ( /, 5X, "WARNING: kgrid is set to zero in the wavefunction file.", &
             /, 14X, "The resulting file will only be usable as the fine grid in inteqp.", / )

END SUBROUTINE write_wfng

!-------------------------------------------------------------------------------

SUBROUTINE real_wfng ( ik, ngkdist_l, nb, ns, energy, wfng_dist )

  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode
  USE mp, ONLY : mp_sum
  USE mp_world, ONLY : world_comm

  IMPLICIT NONE

  integer, intent (in) :: ik, ngkdist_l, nb, ns
  real (DP), intent (in) :: energy ( :, : ) ! ( nb, ns )
  complex (DP), intent (inout) :: wfng_dist ( :, :, : ) ! ( ngkdist_l, nb, ns )

  real (DP), PARAMETER :: eps2 = 1.0D-2
  real (DP), PARAMETER :: eps5 = 1.0D-5
  real (DP), PARAMETER :: eps6 = 1.0D-6

  character :: tmpstr*80
  integer :: i, j, k, is, ib, jb, ig, inum, deg, mdeg, inc
  integer :: dimension_span, reduced_span, ierr
  real (DP) :: x
  integer, allocatable :: imap ( :, : )
  integer, allocatable :: inums ( : )
  integer, allocatable :: inull ( : )
  integer, allocatable :: null_map ( :, : )
  real (DP), allocatable :: psi ( :, : )
  real (DP), allocatable :: phi ( :, : )
  real (DP), allocatable :: vec ( : )
  complex (DP), allocatable :: wfc ( : )

  mdeg = 1
  DO is = 1, ns
    DO ib = 1, nb - 1
      deg = 1
      DO jb = ib + 1, nb
        IF ( abs ( energy ( ib, is ) - energy ( jb, is ) ) &
          .LT. eps5 * dble ( jb - ib + 1 ) ) deg = deg + 1
      ENDDO
      IF ( deg .GT. mdeg ) mdeg = deg
    ENDDO
  ENDDO
  mdeg = mdeg * 2

  ALLOCATE ( imap ( nb, ns ) )
  ALLOCATE ( inums ( ns ) )
  ALLOCATE ( inull ( nb ) )
  ALLOCATE ( null_map ( mdeg, nb ) )

  DO is = 1, ns
    inum = 1
    DO ib = 1, nb
      IF ( ib .EQ. nb ) THEN
        imap ( inum, is ) = ib
        inum = inum + 1
      ELSEIF ( abs ( energy ( ib, is ) - &
        energy ( ib + 1, is ) ) .GT. eps5 ) THEN
        imap ( inum, is ) = ib
        inum = inum + 1
      ENDIF
    ENDDO
    inum = inum - 1
    inums ( is ) = inum
  ENDDO

  ALLOCATE ( wfc ( ngkdist_l ) )
  ALLOCATE ( psi ( ngkdist_l, mdeg ) )
  ALLOCATE ( phi ( ngkdist_l, mdeg ) )
  ALLOCATE ( vec ( ngkdist_l ) )

  DO is = 1, ns
    inc = 1
    inum = inums ( is )
    DO i = 1, inum
      inull ( i ) = 1
      DO ib = inc, imap ( i, is )
        DO ig = 1, ngkdist_l
          wfc ( ig ) = wfng_dist ( ig, ib, is )
        ENDDO
        x = 0.0D0
        DO ig = 1, ngkdist_l
          x = x + dble ( wfc ( ig ) ) **2
        ENDDO
        CALL mp_sum ( x, world_comm )
        IF ( x .LT. eps2 ) null_map ( inull ( i ), i ) = 0
        IF ( x .GT. eps2 ) null_map ( inull ( i ), i ) = 1
        inull ( i ) = inull ( i ) + 1
        x = 0.0D0
        DO ig = 1, ngkdist_l
          x = x + aimag ( wfc ( ig ) ) **2
        ENDDO
        CALL mp_sum ( x, world_comm )
        IF ( x .LT. eps2 ) null_map ( inull ( i ), i ) = 0
        IF ( x .GT. eps2 ) null_map ( inull ( i ), i ) = 1
        inull ( i ) = inull ( i ) + 1
      ENDDO
      inull ( i ) = inull ( i ) - 1
      inc = imap ( i, is ) + 1
    ENDDO
    inc = 1
    ib = 1
    DO i = 1, inum
      k = 1
      DO j = 1, 2 * ( imap ( i, is ) - inc ) + 1, 2
        IF ( null_map ( j, i ) .EQ. 1 .OR. &
          null_map ( j + 1, i ) .EQ. 1 ) THEN
          DO ig = 1, ngkdist_l
            wfc ( ig ) = wfng_dist ( ig, ib, is )
          ENDDO
          IF ( null_map ( j, i ) .EQ. 1 ) THEN
            DO ig = 1, ngkdist_l
              phi ( ig, k ) = dble ( wfc ( ig ) )
            ENDDO
            k = k + 1
          ENDIF
          IF ( null_map ( j + 1, i ) .EQ. 1 ) THEN
            DO ig = 1, ngkdist_l
              phi ( ig, k ) = aimag ( wfc ( ig ) )
            ENDDO
            k = k + 1
          ENDIF
          ib = ib + 1
        ENDIF
      ENDDO
      dimension_span = k - 1
      IF ( dimension_span .EQ. 0 ) THEN
        ierr = 201
        WRITE ( tmpstr, 201 ) ik, is, inc
        CALL errore ( 'real_wfng', tmpstr, ierr )
      ENDIF
      DO j = 1, dimension_span
        x = 0.0D0
        DO ig = 1, ngkdist_l
          x = x + phi ( ig, j ) **2
        ENDDO
        CALL mp_sum ( x, world_comm )
        x = sqrt ( x )
        DO ig = 1, ngkdist_l
          phi ( ig, j ) = phi ( ig, j ) / x
        ENDDO
      ENDDO
!
! the Gram-Schmidt process begins
!
      reduced_span = 1
      DO ig = 1, ngkdist_l
        psi ( ig, 1 ) = phi ( ig, 1 )
      ENDDO
      DO j = 1, dimension_span - 1
        DO ig = 1, ngkdist_l
          vec ( ig ) = phi ( ig, j + 1 )
        ENDDO
        DO k = 1, reduced_span
          x = 0.0D0
          DO ig = 1, ngkdist_l
            x = x + phi ( ig, j + 1 ) * psi ( ig, k )
          ENDDO
          CALL mp_sum ( x, world_comm )
          DO ig = 1, ngkdist_l
            vec ( ig ) = vec ( ig ) - psi ( ig, k ) * x
          ENDDO
        ENDDO
        x = 0.0D0
        DO ig = 1, ngkdist_l
          x = x + vec ( ig ) **2
        ENDDO
        CALL mp_sum ( x, world_comm )
        x = sqrt ( x )
        IF ( x .GT. eps6 ) THEN
          reduced_span = reduced_span + 1
          DO ig = 1, ngkdist_l
            psi ( ig, reduced_span ) = vec ( ig ) / x
          ENDDO
        ENDIF
      ENDDO
!
! the Gram-Schmidt process ends
!
      IF ( reduced_span .LT. imap ( i, is ) - inc + 1 ) THEN
        ierr = 202
        WRITE ( tmpstr, 202 ) ik, is, inc
        CALL errore ( 'real_wfng', tmpstr, ierr )
      ENDIF
      DO ib = inc, imap ( i, is )
        DO ig = 1, ngkdist_l
          wfng_dist ( ig, ib, is ) = &
          CMPLX ( psi ( ig, ib - inc + 1 ), 0.0D0, KIND=dp )
        ENDDO
      ENDDO
      inc = imap ( i, is ) + 1
    ENDDO
  ENDDO

  DEALLOCATE ( vec )
  DEALLOCATE ( phi )
  DEALLOCATE ( psi )
  DEALLOCATE ( wfc )
  DEALLOCATE ( null_map )
  DEALLOCATE ( inull )
  DEALLOCATE ( inums )
  DEALLOCATE ( imap )

  RETURN

201 FORMAT("failed Gram-Schmidt dimension span for kpoint =",i6," spin =",i2," band =",i6)
202 FORMAT("failed Gram-Schmidt reduced span for kpoint =",i6," spin =",i2," band =",i6)

END SUBROUTINE real_wfng

!-------------------------------------------------------------------------------

SUBROUTINE write_rhog ( output_file_name, real_or_complex, symm_type, &
  rhog_nvmin, rhog_nvmax )

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  USE constants, ONLY : pi, tpi, eps6
  USE fft_base, ONLY : dfftp
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, ecutrho
  USE io_global, ONLY : ionode
  USE ions_base, ONLY : nat, atm, ityp, tau 
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rhoz_or_updw
  USE symm_base, ONLY : s, ftau, nsym
  USE matrix_inversion

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  integer, intent (in) :: real_or_complex
  character ( len = 9 ), intent (in) :: symm_type
  integer, intent (in) :: rhog_nvmin
  integer, intent (in) :: rhog_nvmax

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  integer :: unit, id, is, ig, i, j, k, ierr
  integer :: nd, ns, ng_l, ng_g
  integer :: ntran, cell_symmetry, nrecord
  real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: g_g ( :, : )
  complex (DP), allocatable :: rhog_g ( :, : )

  INTEGER, EXTERNAL :: atomic_number

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  IF ( real_or_complex .EQ. 1 ) THEN
    WRITE ( stitle, '("RHO-Real",24X)' )
  ELSE
    WRITE ( stitle, '("RHO-Complex",21X)' )
  ENDIF

  unit = 4
  nrecord = 1
  nd = 3

  ns = nspin
  ng_l = ngm
  ng_g = ngm_g

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_rhog', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2 )
    t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( dfftp%nr1 )
    t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( dfftp%nr2 )
    t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( dfftp%nr3 )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  IF ( rhog_nvmin .NE. 0 .AND. rhog_nvmax .NE. 0 ) &
    CALL calc_rhog ( rhog_nvmin, rhog_nvmax )
    !
    IF ( nspin==2 ) CALL rhoz_or_updw(rho, 'r_and_g', 'updw_rhoz')
    !
  ALLOCATE ( g_g ( nd, ng_g ) )
  ALLOCATE ( rhog_g ( ng_g, ns ) )

  DO ig = 1, ng_g
    DO id = 1, nd
      g_g ( id, ig ) = 0
    ENDDO
  ENDDO
  DO is = 1, ns
    DO ig = 1, ng_g
      rhog_g ( ig, is ) = ( 0.0D0, 0.0D0 )
    ENDDO
  ENDDO

  DO ig = 1, ng_l
    g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  DO is = 1, ns
    DO ig = 1, ng_l
      rhog_g ( ig_l2g ( ig ), is ) = rho%of_g ( ig, is )
    ENDDO
  ENDDO

  CALL mp_sum ( g_g, intra_bgrp_comm )
  CALL mp_sum ( rhog_g, intra_bgrp_comm )

  DO is = 1, ns
    DO ig = 1, ng_g
      rhog_g ( ig, is ) = rhog_g ( ig, is ) * CMPLX ( omega, 0.0D0, KIND=dp )
    ENDDO
  ENDDO

  IF ( ionode ) THEN
    OPEN ( unit = unit, file = TRIM ( output_file_name ), &
      form = 'unformatted', status = 'replace' )
    WRITE ( unit ) stitle, sdate, stime
    WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho
    WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    IF ( real_or_complex .EQ. 1 ) THEN
      WRITE ( unit ) ( ( dble ( rhog_g ( ig, is ) ), &
        ig = 1, ng_g ), is = 1, ns )
    ELSE
      WRITE ( unit ) ( ( rhog_g ( ig, is ), &
        ig = 1, ng_g ), is = 1, ns )
    ENDIF
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( rhog_g )
  DEALLOCATE ( g_g )

  RETURN

END SUBROUTINE write_rhog

!-------------------------------------------------------------------------------

SUBROUTINE calc_rhog (rhog_nvmin, rhog_nvmax)

! calc_rhog    Originally By Brad D. Malone    Last Modified (night before his thesis defense)
! Computes charge density by summing over a subset of occupied bands

  USE cell_base, ONLY : omega, tpiba2
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE gvect, ONLY : ngm, g
  USE io_files, ONLY : nwordwfc, iunwfc
  USE klist, ONLY : xk, nkstot, ngk, nks, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : nbgrp, inter_bgrp_comm
  USE noncollin_module, ONLY : nspin_mag
  USE scf, ONLY : rho
  USE symme, ONLY : sym_rho, sym_rho_init
  USE wavefunctions, ONLY : evc, psic
  USE wvfct, ONLY : wg

  IMPLICIT NONE

  integer, intent (in) :: rhog_nvmin
  integer, intent (in) :: rhog_nvmax
  integer, external :: global_kpoint_index
  integer :: npw, ik, is, ib, ig, ir, iks, ike

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  CALL weights ()

  rho%of_r (:, :) = 0.0D0

  ! take psi to R-space, compute rho in R-space
  DO ik = iks, ike
    is = isk (ik)
    npw = ngk ( ik - iks + 1 )
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    DO ib = rhog_nvmin, rhog_nvmax
      psic (:) = (0.0D0, 0.0D0)
      DO ig = 1, npw
        psic (dfftp%nl (igk_k (ig, ik-iks+1))) = evc (ig, ib)
      ENDDO
      CALL invfft ('Rho', psic, dfftp)
      DO ir = 1, dfftp%nnr
        rho%of_r (ir, is) = rho%of_r (ir, is) + wg (ib, ik) / omega &
          * (dble (psic (ir)) **2 + aimag (psic (ir)) **2)
      ENDDO
    ENDDO
  ENDDO
  CALL mp_sum (rho%of_r, inter_pool_comm)
  CALL mp_sum (rho%of_r, inter_bgrp_comm)
  rho%of_r = rho%of_r / nbgrp

  ! take rho to G-space
  DO is = 1, nspin
    psic (:) = (0.0D0, 0.0D0)
    psic (:) = rho%of_r (:, is)
    CALL fwfft ('Rho', psic, dfftp)
    rho%of_g (:, is) = psic (dfftp%nl (:))
  ENDDO

  ! symmetrize rho (didn`t make a difference)
  CALL sym_rho_init (.False.)
  CALL sym_rho (nspin_mag, rho%of_g)

  RETURN

END SUBROUTINE calc_rhog

!-------------------------------------------------------------------------------

SUBROUTINE write_vxcg ( output_file_name, real_or_complex, symm_type, &
  vxc_zero_rho_core )

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  USE constants, ONLY : pi, tpi, eps6
  USE ener, ONLY : etxc, vtxc
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, ecutrho
  USE io_global, ONLY : ionode
  USE ions_base, ONLY : nat, atm, ityp, tau 
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE symm_base, ONLY : s, ftau, nsym
  USE wavefunctions, ONLY : psic
  USE matrix_inversion

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  integer, intent (in) :: real_or_complex
  character ( len = 9 ), intent (in) :: symm_type
  logical, intent (in) :: vxc_zero_rho_core

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  integer :: unit, id, is, ir, ig, i, j, k, ierr
  integer :: nd, ns, nr, ng_l, ng_g
  integer :: ntran, cell_symmetry, nrecord
  real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: g_g ( :, : )
  real (DP), allocatable :: vxcr_g ( :, : )
  complex (DP), allocatable :: vxcg_g ( :, : )

  INTEGER, EXTERNAL :: atomic_number

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  IF ( real_or_complex .EQ. 1 ) THEN
    WRITE ( stitle, '("VXC-Real",24X)' )
  ELSE
    WRITE ( stitle, '("VXC-Complex",21X)' )
  ENDIF

  unit = 4
  nrecord = 1
  nd = 3

  ns = nspin
  nr = dfftp%nnr
  ng_l = ngm
  ng_g = ngm_g

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_vxcg', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2 )
    t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( dfftp%nr1 )
    t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( dfftp%nr2 )
    t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( dfftp%nr3 )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  ALLOCATE ( g_g ( nd, ng_g ) )
  ALLOCATE ( vxcr_g ( nr, ns ) )
  ALLOCATE ( vxcg_g ( ng_g, ns ) )

  DO ig = 1, ng_g
    DO id = 1, nd
      g_g ( id, ig ) = 0
    ENDDO
  ENDDO
  DO is = 1, ns
    DO ig = 1, ng_g
      vxcg_g ( ig, is ) = ( 0.0D0, 0.0D0 )
    ENDDO
  ENDDO

  DO ig = 1, ng_l
    g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  vxcr_g ( :, : ) = 0.0D0
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  !
  CALL v_xc ( rho, rho_core, rhog_core, etxc, vtxc, vxcr_g )
  !
  DO is = 1, ns
    DO ir = 1, nr
      psic ( ir ) = CMPLX ( vxcr_g ( ir, is ), 0.0D0, KIND=dp )
    ENDDO
    CALL fwfft ( 'Rho', psic, dfftp )
    DO ig = 1, ng_l
      vxcg_g ( ig_l2g ( ig ), is ) = psic ( dfftp%nl ( ig ) )
    ENDDO
  ENDDO

  CALL mp_sum ( g_g, intra_bgrp_comm )
  CALL mp_sum ( vxcg_g, intra_bgrp_comm )

  IF ( ionode ) THEN
    OPEN ( unit = unit, file = TRIM ( output_file_name ), &
      form = 'unformatted', status = 'replace' )
    WRITE ( unit ) stitle, sdate, stime
    WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho
    WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    IF ( real_or_complex .EQ. 1 ) THEN
      WRITE ( unit ) ( ( dble ( vxcg_g ( ig, is ) ), &
        ig = 1, ng_g ), is = 1, ns )
    ELSE
      WRITE ( unit ) ( ( vxcg_g ( ig, is ), &
        ig = 1, ng_g ), is = 1, ns )
    ENDIF
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( vxcg_g )
  DEALLOCATE ( vxcr_g )
  DEALLOCATE ( g_g )

  RETURN

END SUBROUTINE write_vxcg

!-------------------------------------------------------------------------------

SUBROUTINE write_vxc0 ( output_file_name, vxc_zero_rho_core )

  USE constants, ONLY : RYTOEV
  USE ener, ONLY : etxc, vtxc
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE gvect, ONLY : ngm, mill
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE wavefunctions, ONLY : psic

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  logical, intent (in) :: vxc_zero_rho_core

  integer :: unit
  integer :: is, ir, ig
  integer :: nd, ns, nr, ng_l
  real (DP), allocatable :: vxcr_g ( :, : )
  complex (DP), allocatable :: vxc0_g ( : )

  unit = 4
  nd = 3

  ns = nspin
  nr = dfftp%nnr
  ng_l = ngm

  ALLOCATE ( vxcr_g ( nr, ns ) )
  ALLOCATE ( vxc0_g ( ns ) )

  DO is = 1, ns
    vxc0_g ( is ) = ( 0.0D0, 0.0D0 )
  ENDDO

  vxcr_g ( :, : ) = 0.0D0
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  !
  CALL v_xc ( rho, rho_core, rhog_core, etxc, vtxc, vxcr_g )
  !
  DO is = 1, ns
    DO ir = 1, nr
      psic ( ir ) = CMPLX ( vxcr_g ( ir, is ), 0.0D0, KIND=dp )
    ENDDO
    CALL fwfft ( 'Rho', psic, dfftp )
    DO ig = 1, ng_l
      IF ( mill ( 1, ig ) .EQ. 0 .AND. mill ( 2, ig ) .EQ. 0 .AND. &
        mill ( 3, ig ) .EQ. 0 ) vxc0_g ( is ) = psic ( dfftp%nl ( ig ) )
    ENDDO
  ENDDO

  CALL mp_sum ( vxc0_g, intra_bgrp_comm )

  DO is = 1, ns
    vxc0_g ( is ) = vxc0_g ( is ) * CMPLX ( RYTOEV, 0.0D0, KIND=dp )
  ENDDO

  IF ( ionode ) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    WRITE ( unit, 101 )
    DO is = 1, ns
      WRITE ( unit, 102 ) is, vxc0_g ( is )
    ENDDO
    WRITE ( unit, 103 )
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  DEALLOCATE ( vxcr_g )
  DEALLOCATE ( vxc0_g )

  RETURN

101 FORMAT ( /, 5X, "--------------------------------------------", &
             /, 5X, "spin    Re Vxc(G=0) (eV)    Im Vxc(G=0) (eV)", &
             /, 5X, "--------------------------------------------" )
102 FORMAT ( 5X, I1, 3X, 2F20.15 )
103 FORMAT ( 5X, "--------------------------------------------", / )

END SUBROUTINE write_vxc0

!-------------------------------------------------------------------------------

SUBROUTINE write_vxc_r (output_file_name, diag_nmin, diag_nmax, &
  offdiag_nmin, offdiag_nmax, vxc_zero_rho_core)

  USE kinds, ONLY : DP
  USE constants, ONLY : rytoev
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : invfft
  USE gvect, ONLY : ngm, g
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global,        ONLY : ionode, ionode_id, stdout    !FZ: test if need this when output to pp_out
  !USE io_global, ONLY : ionode       !FZ: comment
  USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE wavefunctions, ONLY : evc, psic
  USE wvfct, ONLY : nbnd

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax
  integer, intent (inout) :: offdiag_nmin
  integer, intent (inout) :: offdiag_nmax
  logical, intent (in) :: vxc_zero_rho_core

  integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2
  integer, external :: global_kpoint_index
  real (DP) :: dummyr
  complex (DP) :: dummyc
  real (DP), allocatable :: mtxeld (:, :)
  complex (DP), allocatable :: mtxelo (:, :, :)
  real (DP), allocatable :: vxcr (:, :)
  complex (DP), allocatable :: psic2 (:)

  if(diag_nmin > diag_nmax) then
    call errore ( 'write_vxc_r', 'diag_nmin > diag_nmax', diag_nmin )
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    diag_nmax = nbnd
  ENDIF
  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

  if(offdiag_nmin > offdiag_nmax) then
    call errore ( 'write_vxc_r', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
  endif
  IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
  IF (offdiag_nmax .GT. nbnd)  then
    write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
    offdiag_nmax = nbnd
  ENDIF
  noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

  IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

  unit = 4

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = 0.0D0
  ENDIF
  IF (noffdiag .GT. 0) THEN
    ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
    mtxelo (:, :, :) = (0.0D0, 0.0D0)
  ENDIF

  ALLOCATE (vxcr (dfftp%nnr, nspin))
  IF (noffdiag .GT. 0) ALLOCATE (psic2 (dfftp%nnr))

  vxcr (:, :) = 0.0D0
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  !
  CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxcr)
  !
  IF ( ionode ) WRITE (stdout, '(" Test that stdout write_vxc_r is called")')    !FZ: test
  IF ( ionode ) WRITE (6, '(" Test that WRITE(6, ...) write_vxc_r is called")')  !FZ: test
  !IF (ionode) THEN                                          !FZ: test
  !  OPEN (unit = unit, file = TRIM (output_file_name) // 'test_pw2bgw_output_file', &  !FZ: test
  !    form = 'formatted', status = 'replace')               !FZ: test
  DO ik = iks, ike
    npw = ngk ( ik - iks + 1 )
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    !WRITE (unit, 102) ndiag, noffdiag, ndiag, noffdiag           !FZ test
    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (dfftp%nl (igk_k (ig,ik-iks+1))) = evc (ig, ib)
        ENDDO
        CALL invfft ('Rho', psic, dfftp)
        dummyr = 0.0D0
        DO ir = 1, dfftp%nnr
          dummyr = dummyr + vxcr (ir, isk (ik)) &
            * (dble (psic (ir)) **2 + aimag (psic (ir)) **2)
        ENDDO
        dummyr = dummyr * rytoev / dble (dfftp%nr1x * dfftp%nr2x * dfftp%nr3x)
        CALL mp_sum (dummyr, intra_bgrp_comm)
        mtxeld (ib - diag_nmin + 1, ik) = dummyr
      ENDDO
      !WRITE (unit, 102) is, ib, dummyr, dble (dfftp%nr1x * dfftp%nr2x * dfftp%nr3x)    !FZ  test
    ENDIF
    IF (noffdiag .GT. 0) THEN
      DO ib = offdiag_nmin, offdiag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (dfftp%nl (igk_k (ig,ik-iks+1))) = evc (ig, ib)
        ENDDO
        CALL invfft ('Rho', psic, dfftp)
        DO ib2 = offdiag_nmin, offdiag_nmax
          psic2 (:) = (0.0D0, 0.0D0)
          DO ig = 1, npw
            psic2 (dfftp%nl (igk_k (ig,ik-iks+1))) = evc (ig, ib2)
          ENDDO
          CALL invfft ('Rho', psic2, dfftp)
          dummyc = (0.0D0, 0.0D0)
          DO ir = 1, dfftp%nnr
            dummyc = dummyc + CMPLX (vxcr (ir, isk (ik)), 0.0D0, KIND=dp) &
              * conjg (psic2 (ir)) * psic (ir)
          ENDDO
          dummyc = dummyc &
               * CMPLX (rytoev / dble (dfftp%nr1x * dfftp%nr2x * dfftp%nr3x), &
                        0.0D0, KIND=dp)
          CALL mp_sum (dummyc, intra_bgrp_comm)
          mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin &
            + 1, ik) = dummyc
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  !  CLOSE (unit = unit, status = 'keep')      !FZ:  test
  !ENDIF                                       !FZ:  test

  DEALLOCATE (vxcr)
  IF (noffdiag .GT. 0) DEALLOCATE (psic2)

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
  IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
      WRITE (unit, '(" Test that write_vxc_r is called")')    !FZ:  test
      WRITE (unit, 101) xk(:, 1), nspin * ndiag, nspin * noffdiag    !FZ:  test
    DO ik = 1, nkstot / nspin
      WRITE (unit, 101) xk(:, ik), nspin * ndiag, &
        nspin * noffdiag **2
      DO is = 1, nspin
        IF (ndiag .GT. 0) THEN
          DO ib = diag_nmin, diag_nmax
            WRITE (unit, 102) is, ib, mtxeld &
              (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin), &
              0.0D0
          ENDDO
        ENDIF
        IF (noffdiag .GT. 0) THEN
          DO ib = offdiag_nmin, offdiag_nmax
            DO ib2 = offdiag_nmin, offdiag_nmax
              WRITE (unit, 103) is, ib2, ib, mtxelo &
                (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                ik + (is - 1) * nkstot / nspin)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
  IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

  RETURN

  101 FORMAT (3F13.9, 2I8)
  102 FORMAT (2I8, 2F15.9)
  103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_vxc_r

!-------------------------------------------------------------------------------

SUBROUTINE write_v_h_g (output_file_name, diag_nmin, diag_nmax, &     !FZ: all function is for metaGGA
  offdiag_nmin, offdiag_nmax)                      !FZ: output Hartree energy for each bands and kpts

  USE constants, ONLY : rytoev
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc, ehart !FZ:  added ehart
  USE exx, ONLY : vexx
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  !USE funct, ONLY : exx_is_active      !FZ: commented
  !USE funct, ONLY : exx_is_active, dft_is_meta, get_meta     !FZ: Added for metaGGA commented in v_h
  USE gvect, ONLY : ngm, g
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho !, rho_core, rhog_core  FZ: for vhartree
  USE wavefunctions, ONLY : evc, psic
  USE wvfct, ONLY : npwx, nbnd

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax
  integer, intent (inout) :: offdiag_nmin
  integer, intent (inout) :: offdiag_nmax

  integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2, ikk
  integer, external :: global_kpoint_index
  complex (DP) :: dummy
  complex (DP), allocatable :: mtxeld (:, :)
  complex (DP), allocatable :: mtxelo (:, :, :)
  real (DP), allocatable :: vxcr (:, :)
  real (DP), allocatable :: v_h (:, :)       !FZ: for metaGGA , stores hartree potential
  real (DP), allocatable :: kedtaur (:, :)   !FZ: for metaGGA   kedtau: local K energy density,  kedtaur (in realspace)
  complex (DP), allocatable :: psic2 (:)
  complex (DP), allocatable :: hpsi (:)
  REAL(DP)  ::  charge       !FZ: for metaGGA, stores total hartree energy, charge

  if(diag_nmin > diag_nmax) then
    call errore ( 'write_v_h_g', 'diag_nmin > diag_nmax', diag_nmin )  !FZ:
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    diag_nmax = nbnd
  ENDIF
  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

  if(offdiag_nmin > offdiag_nmax) then
    call errore ( 'write_v_h_g', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )  !FZ:
  endif
  IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
  IF (offdiag_nmax .GT. nbnd)  then
    write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
    offdiag_nmax = nbnd
  ENDIF
  noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

  IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

  unit = 4

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = (0.0D0, 0.0D0)
  ENDIF
  IF (noffdiag .GT. 0) THEN
    ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
    mtxelo (:, :, :) = (0.0D0, 0.0D0)
  ENDIF

  ALLOCATE (vxcr (dfftp%nnr, nspin))
  ALLOCATE (v_h (dfftp%nnr, nspin))     !FZ: for meta GGA
  IF (noffdiag .GT. 0) ALLOCATE (psic2 (dfftp%nnr))
  ALLOCATE (hpsi (dfftp%nnr))
  ALLOCATE (kedtaur (dfftp%nnr, nspin))   !FZ: for meta GGA
  vxcr (:, :) = 0.0D0
  v_h (:, :) = 0.0D0        !FZ: for metaGGA
  kedtaur (:, :) = 0.0D0    !FZ: for metaGGA
  ehart  = 0.D0
  charge = 0.D0
  CALL v_h_only (rho%of_g(:,1), ehart, charge, v_h)   !FZ: for metaGGA   important check!rho%of_g(:,1)
  
  IF (ionode) THEN                                          !FZ: test
    OPEN (unit = unit, file = 'test_v_h_output', &          !FZ: test
      form = 'formatted', status = 'replace')               !FZ: test
      WRITE (unit, *) "v_h = ", v_h  !FZ:  test
    CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  ENDIF                                          !FZ:  test

  DO ik = iks, ike
    ikk = ik - iks + 1
    npw = ngk ( ik - iks + 1 )
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (dfftp%nl (igk_k(ig,ikk))) = evc (ig, ib)      !FZ:   evc stores all eigenfunctions in G space, use davcio function to extract ik - iks + 1 th kpoint 's all the eigen functions for all bands, 
   !FZ:   e.g.  silicon 10 kpts, 22 plane waves, so ig can be 1-22, igk_k(ig, ikk) reorders # of plane waves for each k point ikk
   !FZ:  igk_k is a 22*10 matrix
        ENDDO
        CALL invfft ('Rho', psic, dfftp)
        DO ir = 1, dfftp%nnr
          psic (ir) = psic (ir) * v_h (ir, isk (ik))      !FZ:   v_hartree | psi>
        ENDDO
        CALL fwfft ('Rho', psic, dfftp)
        hpsi (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          hpsi (ig) = psic (dfftp%nl (igk_k(ig,ikk)))
        ENDDO
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (ig) = evc (ig, ib)
        ENDDO
        !IF (exx_is_active ()) &
        !   CALL vexx (npwx, npw, 1, psic, hpsi)
        dummy = (0.0D0, 0.0D0)
        DO ig = 1, npw
          dummy = dummy + conjg (psic (ig)) * hpsi (ig)
        ENDDO
        IF (ionode) THEN                                          !FZ: test
          OPEN (unit = unit, file = 'test_v_h_output', &  !FZ: test
            form = 'formatted', status = 'replace')               !FZ: test
            WRITE (unit, *) "is = ",is," ib = ", ib, " dummy = ",dummy      !FZ test
          CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
        ENDIF                                       !FZ:  test
        dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
        CALL mp_sum (dummy, intra_bgrp_comm)
        mtxeld (ib - diag_nmin + 1, ik) = dummy
      ENDDO
    ENDIF
    IF (noffdiag .GT. 0) THEN
      DO ib = offdiag_nmin, offdiag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (dfftp%nl (igk_k(ig,ikk))) = evc (ig, ib)
        ENDDO
        CALL invfft ('Rho', psic, dfftp)
        DO ir = 1, dfftp%nnr
          psic (ir) = psic (ir) * v_h (ir, isk (ik))     !FZ: for metaGGA
        ENDDO
        CALL fwfft ('Rho', psic, dfftp)
        hpsi (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          hpsi (ig) = psic (dfftp%nl (igk_k (ig,ikk)))
        ENDDO
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (ig) = evc (ig, ib)
        ENDDO
        !IF (exx_is_active ()) &
        !   CALL vexx (npwx, npw, 1, psic, hpsi)
        DO ib2 = offdiag_nmin, offdiag_nmax
          psic2 (:) = (0.0D0, 0.0D0)
          DO ig = 1, npw
            psic2 (ig) = evc (ig, ib2)
          ENDDO
          dummy = (0.0D0, 0.0D0)
          DO ig = 1, npw
            dummy = dummy + conjg (psic2 (ig)) * hpsi (ig)
          ENDDO
          dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
          CALL mp_sum (dummy, intra_bgrp_comm)
          mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin &
            + 1, ik) = dummy
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  IF (ionode) THEN                                          !FZ: test
    OPEN (unit = unit, file = 'test_v_h_mtxeld', &  !FZ: test
      form = 'formatted', status = 'replace')               !FZ: test
      WRITE (unit, *) "mtxeld = ", mtxeld  !FZ:  test
      !WRITE (unit, *) "exx_is_active = ",exx_is_active ()  !FZ:  test
    CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  ENDIF                                       !FZ:  test

  DEALLOCATE (vxcr)
  DEALLOCATE (v_h)          !FZ: for metaGGA calculate hartree
  DEALLOCATE (kedtaur)      !FZ: for metaGGA
  IF (noffdiag .GT. 0) DEALLOCATE (psic2)
  DEALLOCATE (hpsi)

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
  IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    DO ik = 1, nkstot / nspin
      WRITE (unit, 101) xk(:, ik), nspin * ndiag, &
        nspin * noffdiag **2
      DO is = 1, nspin
        IF (ndiag .GT. 0) THEN
          DO ib = diag_nmin, diag_nmax
            WRITE (unit, 102) is, ib, mtxeld &
              (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin)
          ENDDO
        ENDIF
        IF (noffdiag .GT. 0) THEN
          DO ib = offdiag_nmin, offdiag_nmax
            DO ib2 = offdiag_nmin, offdiag_nmax
              WRITE (unit, 103) is, ib2, ib, mtxelo &
                (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                ik + (is - 1) * nkstot / nspin)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
  IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

  RETURN

  101 FORMAT (3F13.9, 2I8)
  102 FORMAT (2I8, 2F15.9)
  103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_v_h_g

!-------------------------------------------------------------------------------

SUBROUTINE write_vltot (output_file_name, diag_nmin, diag_nmax, &     !FZ: all function is for metaGGA
  offdiag_nmin, offdiag_nmax)                      !FZ: output local (ionic) potential energy for each bands and kpts

  USE constants, ONLY : rytoev
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc, ehart !FZ:  added ehart
  USE exx, ONLY : vexx
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  !USE funct, ONLY : exx_is_active      !FZ: commented
  !USE funct, ONLY : exx_is_active, dft_is_meta, get_meta     !FZ: Added for metaGGA commented in v_h
  USE gvect, ONLY : ngm, g ,gstart  !FZ: added gstart 
  USE gvecs, ONLY : doublegrid   !FZ: added for calculation of vrs
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, vltot, vrs, v, kedtau !, rho_core, rhog_core  FZ: vltot for total local potential (vltot(:,:) is in real space)
  USE wavefunctions, ONLY : evc, psic
  USE wvfct, ONLY : npwx, nbnd, et

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax
  integer, intent (inout) :: offdiag_nmin
  integer, intent (inout) :: offdiag_nmax

  integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2, ikk
  integer, external :: global_kpoint_index
  complex (DP) :: dummy
  complex (DP), allocatable :: mtxeld (:, :)
  complex (DP), allocatable :: mtxelo (:, :, :)
  !real (DP), allocatable :: vxcr (:, :)
  !real (DP), allocatable :: v_h (:, :)       !FZ: for metaGGA , stores hartree potential
  !real (DP), allocatable :: kedtaur (:, :)   !FZ: for metaGGA   kedtau: local K energy density,  kedtaur (in realspace)
  complex (DP), allocatable :: psic2 (:)
  complex (DP), allocatable :: hpsi (:)
  !REAL(DP)  ::  charge       !FZ: for metaGGA, stores total hartree energy, charge

  if(diag_nmin > diag_nmax) then
    call errore ( 'write_vltot', 'diag_nmin > diag_nmax', diag_nmin )  !FZ:
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    diag_nmax = nbnd
  ENDIF
  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

  if(offdiag_nmin > offdiag_nmax) then
    call errore ( 'write_vltot', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )  !FZ:
  endif
  IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
  IF (offdiag_nmax .GT. nbnd)  then
    write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
    offdiag_nmax = nbnd
  ENDIF
  noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

  IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

  unit = 4

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = (0.0D0, 0.0D0)
  ENDIF
  IF (noffdiag .GT. 0) THEN
    ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
    mtxelo (:, :, :) = (0.0D0, 0.0D0)
  ENDIF

  !ALLOCATE (vxcr (dfftp%nnr, nspin))
  !ALLOCATE (v_h (dfftp%nnr, nspin))     !FZ: for meta GGA
  IF (noffdiag .GT. 0) ALLOCATE (psic2 (dfftp%nnr))
  ALLOCATE (hpsi (dfftp%nnr))
  !ALLOCATE (kedtaur (dfftp%nnr, nspin))   !FZ: for meta GGA
  !vxcr (:, :) = 0.0D0
  !v_h (:, :) = 0.0D0        !FZ: for metaGGA
  !kedtaur (:, :) = 0.0D0    !FZ: for metaGGA
  !ehart  = 0.D0
  !charge = 0.D0
  !CALL v_h_only (rho%of_g(:,1), ehart, charge, v_h)   !FZ: for metaGGA   important check!rho%of_g(:,1)
  CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, &    !FZ: calculate vrs which is vltot + v_H + vxc
                         nspin, doublegrid )     !FZ: calculate vrs which is vltot + v_H + vxc
  
  IF (ionode) THEN                                          !FZ: test
    OPEN (unit = unit, file = 'test_vltot_output', &          !FZ: test
      form = 'formatted', status = 'replace')               !FZ: test
      WRITE (unit, *) "vltot = ", vltot  !FZ:  test
    CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  ENDIF                                          !FZ:  test
  IF (ionode) THEN                                          !FZ: test
    OPEN (unit = unit, file = 'test_vrs_output', &          !FZ: test
      form = 'formatted', status = 'replace')               !FZ: test
      WRITE (unit, *) "vrs = ", vrs  !FZ:  test
    CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  ENDIF                                          !FZ:  test
  IF (ionode) THEN                                          !FZ: test
    OPEN (unit = unit, file = 'test_et_output', &          !FZ: test
      form = 'formatted', status = 'replace')               !FZ: test
      WRITE (unit, *) "gstart = ", gstart  !FZ:  test
      WRITE (unit, *) "et = ", et * CMPLX (rytoev, 0.0D0, KIND=dp)  !FZ:  test
    CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  ENDIF                                          !FZ:  test

  DO ik = iks, ike
    ikk = ik - iks + 1
    npw = ngk ( ik - iks + 1 )
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (dfftp%nl (igk_k(ig,ikk))) = evc (ig, ib)      !FZ:   evc stores all eigenfunctions in G space, use davcio function to extract ik - iks + 1 th kpoint 's all the eigen functions for all bands, 
   !FZ:   e.g.  silicon 10 kpts, 22 plane waves, so ig can be 1-22, igk_k(ig, ikk) reorders # of plane waves for each k point ikk
   !FZ:  igk_k is a 22*10 matrix
        ENDDO
        CALL invfft ('Rho', psic, dfftp)
        DO ir = 1, dfftp%nnr
          !psic (ir) = psic (ir) * vltot (ir)  !FZ: vltot    !FZ:   v local total | psi>
          psic (ir) = psic (ir) * vrs (ir, isk(ik))  !FZ: vrs    !FZ:   v local total + vr | psi>
          !psic (ir) = psic (ir) * (v%of_r (ir, isk(ik)) + vltot(ir)) !FZ: vsc    !FZ:   v local total + vr | psi>
             !FZ:   vrs (nrxx, nspin), vltot (nrxx), vr (nrxx, nspin)  (v_h has the same structure as vr)  for v_h it is v_h (ir, isk (ik)), for vltot it is vltot(ir)
        ENDDO
        CALL fwfft ('Rho', psic, dfftp)
        hpsi (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          hpsi (ig) = psic (dfftp%nl (igk_k(ig,ikk)))
        ENDDO
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (ig) = evc (ig, ib)
        ENDDO
        !IF (exx_is_active ()) &
        !   CALL vexx (npwx, npw, 1, psic, hpsi)
        dummy = (0.0D0, 0.0D0)
        DO ig = 1, npw
          dummy = dummy + conjg (psic (ig)) * hpsi (ig)
        ENDDO
        IF (ionode) THEN                                          !FZ: test
          OPEN (unit = unit, file = 'test_vltot_output2', &  !FZ: test
            form = 'formatted', status = 'replace')               !FZ: test
            WRITE (unit, *) "is = ",is," ib = ", ib, " dummy = ",dummy      !FZ test
          CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
        ENDIF                                       !FZ:  test
        dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
        CALL mp_sum (dummy, intra_bgrp_comm)
        mtxeld (ib - diag_nmin + 1, ik) = dummy
      ENDDO
    ENDIF
    IF (noffdiag .GT. 0) THEN
      DO ib = offdiag_nmin, offdiag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (dfftp%nl (igk_k(ig,ikk))) = evc (ig, ib)
        ENDDO
        CALL invfft ('Rho', psic, dfftp)
        DO ir = 1, dfftp%nnr
          !psic (ir) = psic (ir) * vltot (ir)  !FZ: vltot    !FZ:   v local total | psi>
          psic (ir) = psic (ir) * vrs (ir, isk(ik))  !FZ: vrs    !FZ:   v local total + vr | psi>
          !psic (ir) = psic (ir) * (v%of_r (ir, isk(ik)) + vltot(ir)) !FZ: vsc    !FZ:   v local total + vr | psi>
        ENDDO
        CALL fwfft ('Rho', psic, dfftp)
        hpsi (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          hpsi (ig) = psic (dfftp%nl (igk_k (ig,ikk)))
        ENDDO
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (ig) = evc (ig, ib)
        ENDDO
        !IF (exx_is_active ()) &
        !   CALL vexx (npwx, npw, 1, psic, hpsi)
        DO ib2 = offdiag_nmin, offdiag_nmax
          psic2 (:) = (0.0D0, 0.0D0)
          DO ig = 1, npw
            psic2 (ig) = evc (ig, ib2)
          ENDDO
          dummy = (0.0D0, 0.0D0)
          DO ig = 1, npw
            dummy = dummy + conjg (psic2 (ig)) * hpsi (ig)
          ENDDO
          dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
          CALL mp_sum (dummy, intra_bgrp_comm)
          mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin &
            + 1, ik) = dummy
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  IF (ionode) THEN                                          !FZ: test
    OPEN (unit = unit, file = 'test_vltot_mtxeld', &  !FZ: test
      form = 'formatted', status = 'replace')               !FZ: test
      WRITE (unit, *) "mtxeld = ", mtxeld  !FZ:  test
      !WRITE (unit, *) "exx_is_active = ",exx_is_active ()  !FZ:  test
    CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  ENDIF                                       !FZ:  test

  !DEALLOCATE (vxcr)
  !DEALLOCATE (v_h)          !FZ: for metaGGA calculate hartree
  !DEALLOCATE (kedtaur)      !FZ: for metaGGA
  IF (noffdiag .GT. 0) DEALLOCATE (psic2)
  DEALLOCATE (hpsi)

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
  IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    DO ik = 1, nkstot / nspin
      WRITE (unit, 101) xk(:, ik), nspin * ndiag, &
        nspin * noffdiag **2
      DO is = 1, nspin
        IF (ndiag .GT. 0) THEN
          DO ib = diag_nmin, diag_nmax
            WRITE (unit, 102) is, ib, mtxeld &
              (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin)
          ENDDO
        ENDIF
        IF (noffdiag .GT. 0) THEN
          DO ib = offdiag_nmin, offdiag_nmax
            DO ib2 = offdiag_nmin, offdiag_nmax
              WRITE (unit, 103) is, ib2, ib, mtxelo &
                (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                ik + (is - 1) * nkstot / nspin)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
  IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

  RETURN

  101 FORMAT (3F13.9, 2I8)
  102 FORMAT (2I8, 2F15.9)
  103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_vltot

!-------------------------------------------------------------------------------

SUBROUTINE write_ekin_g (output_file_name, diag_nmin, diag_nmax)    !FZ: all function is for metaGGA 
!FZ: output kinetic energy for each bands and kpts (only diagonal)

  USE constants, ONLY : rytoev
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc, ehart !FZ:  added ehart
  USE exx, ONLY : vexx
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  !USE funct, ONLY : exx_is_active      !FZ: commented
  !USE funct, ONLY : exx_is_active, dft_is_meta, get_meta     !FZ: Added for metaGGA commented in v_h
  USE gvect, ONLY : ngm, g
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, v_of_0 !, rho_core, rhog_core  FZ: for vhartree   !FZ: v_of_0 is the V(G=0) term, need to add to diagonal term which is g2kin,  because in v_H calculation, we did not calculate V(G=0)
  USE wavefunctions, ONLY : evc, psic
  USE wvfct, ONLY : npwx, nbnd, g2kin   !FZ: added g2kin
  USE g_psi_mod,            ONLY : h_diag, s_diag !FZ: for nonlocal pseudo
  USE noncollin_module, ONLY: noncolin, npol  !FZ: for nonlocal pseudo
  USE uspp,                 ONLY : vkb, nkb, okvan, deeq, qq_at, qq_so, deeq_nc, indv_ijkb0 !FZ: test for nonlocal
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp !FZ: test for nonlocal
  USE wvfct, ONLY: npwx !FZ: test for nonlocal
  USE lsda_mod, ONLY: current_spin !FZ: test for nonlocal
  USE uspp_param, ONLY: upf, nh !FZ: test for nonlocal
  USE spin_orb, ONLY: lspinorb !FZ: test for nonlocal

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax

  integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, ib2, ikk
  integer, external :: global_kpoint_index
  complex (DP) :: dummy
  complex (DP), allocatable :: mtxeld (:, :)
  real (DP), allocatable :: vxcr (:, :)
  !real (DP), allocatable :: v_h (:, :)       !FZ: for metaGGA , stores hartree potential
  !real (DP), allocatable :: kedtaur (:, :)   !FZ: for metaGGA   kedtau: local K energy density,  kedtaur (in realspace)
  complex (DP), allocatable :: psic2 (:)
  complex (DP), allocatable :: hpsi (:)
  !REAL(DP)  ::  charge       !FZ: for metaGGA, stores total hartree energy, charge
  character(len=20) :: Ctemp
  INTEGER :: ierr   !FZ: for nonlocal

  if(diag_nmin > diag_nmax) then
    call errore ( 'write_v_h_g', 'diag_nmin > diag_nmax', diag_nmin )  !FZ:
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    diag_nmax = nbnd
  ENDIF
  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)
  
  IF (ndiag .EQ. 0) RETURN

  unit = 4

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = (0.0D0, 0.0D0)
  ENDIF

  ALLOCATE (vxcr (dfftp%nnr, nspin))
  !ALLOCATE (v_h (dfftp%nnr, nspin))     !FZ: for meta GGA
  ALLOCATE (hpsi (dfftp%nnr))
  !ALLOCATE (kedtaur (dfftp%nnr, nspin))   !FZ: for meta GGA
  vxcr (:, :) = 0.0D0
  !v_h (:, :) = 0.0D0        !FZ: for metaGGA
  !kedtaur (:, :) = 0.0D0    !FZ: for metaGGA
  ehart  = 0.D0
  !charge = 0.D0
  !CALL v_h_only (rho%of_g(:,1), ehart, charge, v_h)   !FZ: for metaGGA   important check!rho%of_g(:,1)
  ALLOCATE( h_diag( npwx, npol ), STAT=ierr )    !FZ: for nonlocal
  IF( ierr /= 0 ) &    !FZ: for nonlocal
     CALL errore( ' diag_bands ', ' cannot allocate h_diag ', ABS(ierr) )     !FZ: for nonlocal
  ALLOCATE( s_diag( npwx, npol ), STAT=ierr )     !FZ: for nonlocal
  IF( ierr /= 0 ) &     !FZ: for nonlocal
     CALL errore( ' diag_bands ', ' cannot allocate s_diag ', ABS(ierr) )     !FZ: for nonlocal
  !ALLOCATE (h_diag (npwx)) !FZ: for nonlocal
  !ALLOCATE (s_diag (npwx)) !FZ: for nonlocal
  h_diag (:,:) = 0.0D0  !FZ: for nonlocal
  s_diag (:,:) = 0.0D0  !FZ: for nonlocal
  
  IF (ionode) THEN                                          !FZ: test
    OPEN (unit = unit, file = 'test_ekin_output', &          !FZ: test
      form = 'formatted', status = 'replace')               !FZ: test
      WRITE (unit, *) "g2kin = ", g2kin  !FZ:  test
    CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  ENDIF                                          !FZ:  test
  IF (ionode) THEN                                          !FZ: test
    OPEN (unit = unit, file = 'test_ekin_h_diag_output', &          !FZ: test
      form = 'formatted', status = 'replace')               !FZ: test
      WRITE (unit, *) "h_diag = ", h_diag  !FZ:  test
      WRITE (unit, *) "s_diag = ", s_diag  !FZ:  test
      WRITE (unit, *) "vkb = ", vkb  !FZ:  test
      WRITE (unit, *) "npol = ", npol  !FZ:  test
    CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  ENDIF                                          !FZ:  test

  DO ik = iks, ike
    ikk = ik - iks + 1
    npw = ngk ( ik - iks + 1 )
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    vkb (:,:) = 0.0D0  !FZ: for nonlocal
    IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb ) !FZ: for nonlocal
    h_diag (:,:) = 0.0D0  !FZ: for nonlocal
    s_diag (:,:) = 0.0D0  !FZ: for nonlocal
    CALL usnldiag( npw, h_diag, s_diag )   !FZ: add nonlocal pseudopotential term to diagonal part of Hamiltonian
    g2kin(:) = 0.0D0  !FZ: important: initialize g2kin
    call g2_kin( ik )    !FZ:  check: if index is ik or ikk
    IF (ionode) THEN              !FZ: test
      write(Ctemp,"(I2)") ik                                          !FZ: test
      OPEN (unit = unit, file = 'test_g2_kin_output' // Trim(AdjustL(Ctemp)) //'.dat' , &  !FZ: test
        form = 'formatted', status = 'replace')               !FZ: test
        WRITE (unit, *) "v_of_0 = ", v_of_0      !FZ test
        WRITE (unit, *) "g2_kin = ", g2kin      !FZ test
        WRITE (unit, *) "vkb = ", vkb  !FZ:  test
      CLOSE (unit = unit, status = 'keep')      !FZ:  test
    ENDIF                                       !FZ:  test
    IF (ionode) THEN              !FZ: test
      write(Ctemp,"(I2)") ik                                          !FZ: test
      OPEN (unit = unit, file = 'test_evc_output' // Trim(AdjustL(Ctemp)) //'.dat' , &  !FZ: test
        form = 'formatted', status = 'replace')               !FZ: test
        WRITE (unit, *) "evc = ",evc      !FZ test
      CLOSE (unit = unit, status = 'keep')      !FZ:  test
    ENDIF                                       !FZ:  test
    IF (abs(g2kin(1)) .lt. 0.0001_dp) THEN     !FZ: test
    IF (ionode) THEN                                          !FZ: test
      OPEN (unit = unit, file = 'test_ekin_h_diag_output2', &          !FZ: test
        form = 'formatted', status = 'replace')               !FZ: test
        WRITE (unit, *) "h_diag = ", h_diag  !FZ:  test
        WRITE (unit, *) "s_diag = ", s_diag  !FZ:  test
      CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
    ENDIF                                          !FZ:  test
    ENDIF                                          !FZ:  test
    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        psic (:) = (0.0D0, 0.0D0)
        hpsi (:) = (0.0D0, 0.0D0)  !FZ: added initialization
        DO ig = 1, npw
          !psic (igk_k(ig,ikk)) = evc (ig, ib)      !FZ:   evc stores all eigenfunctions in G space, use davcio function to extract ik - iks + 1 th kpoint 's all the eigen functions for all bands, 
          psic (ig) = evc (ig, ib)      !FZ:   evc stores all eigenfunctions in G space, use davcio function to extract ik - iks + 1 th kpoint 's all the eigen functions for all bands, 
   !FZ:   e.g.  silicon 10 kpts, 22 plane waves, so ig can be 1-22, igk_k(ig, ikk) reorders # of plane waves for each k point ikk
   !FZ:  igk_k is a 22*10 matrix
        ENDDO
        !CALL invfft ('Rho', psic, dfftp)
        !DO ir = 1, dfftp%nnr
        !  psic (ir) = psic (ir) * v_h (ir, isk (ik))      !FZ:   v_hartree | psi>
        !ENDDO
        !CALL fwfft ('Rho', psic, dfftp)
        !hpsi (:) = (0.0D0, 0.0D0)
        !DO ig = 1, npw
        !  hpsi (ig) = psic (dfftp%nl (igk_k(ig,ikk)))
        !ENDDO
        DO ig = 1, npw                                  !FZ: calculate |k+G|^2 |psi>
          !hpsi (ig) = g2kin(ig) * psic (igk_k(ig,ikk))   !FZ: calculate |k+G|^2 |psi>
          hpsi (ig) = g2kin(ig) * psic (ig)   !FZ: calculate |k+G|^2 |psi>
          !hpsi (ig) = (g2kin(ig) + h_diag(ig, 1)) * psic (igk_k(ig,ikk))   !FZ: calculate |k+G|^2 |psi>
        ENDDO                                           !FZ: calculate |k+G|^2 |psi>
        !hpsi(1) = hpsi(1) + v_of_0 * psic (igk_k(1,ikk))  !FZ: add local potential 
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (ig) = evc (ig, ib)      !FZ: 
          !psic (igk_k(ig,ikk)) = evc (ig, ib)  !FZ: wrong
        ENDDO
        !IF (exx_is_active ()) &
        !   CALL vexx (npwx, npw, 1, psic, hpsi)
        dummy = (0.0D0, 0.0D0)
        DO ig = 1, npw
          dummy = dummy + conjg (psic (ig)) * hpsi (ig)
        ENDDO
        IF (ionode) THEN                                          !FZ: test
          OPEN (unit = unit, file = 'test_ekin_output2', &  !FZ: test
            form = 'formatted', status = 'replace')               !FZ: test
            WRITE (unit, *) "is = ",is," ib = ", ib, " dummy = ",dummy      !FZ test
          CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
        ENDIF                                       !FZ:  test
        dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
        CALL mp_sum (dummy, intra_bgrp_comm)
        mtxeld (ib - diag_nmin + 1, ik) = dummy
      ENDDO
    ENDIF
    
  ENDDO
  IF (ionode) THEN                                          !FZ: test
    OPEN (unit = unit, file = 'test_ekin_mtxeld', &  !FZ: test
      form = 'formatted', status = 'replace')               !FZ: test
      WRITE (unit, *) "mtxeld = ", mtxeld  !FZ:  test
      !WRITE (unit, *) "exx_is_active = ",exx_is_active ()  !FZ:  test
    CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  ENDIF                                       !FZ:  test

  DEALLOCATE (vxcr)
  !DEALLOCATE (v_h)          !FZ: for metaGGA calculate hartree
  !DEALLOCATE (kedtaur)      !FZ: for metaGGA
  DEALLOCATE (hpsi)
  DEALLOCATE (h_diag)  !FZ: for nonlocal
  DEALLOCATE (s_diag)  !FZ: for nonlocal

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    DO ik = 1, nkstot / nspin
      WRITE (unit, 101) xk(:, ik), nspin * ndiag!, &   !FZ: commented
        !nspin * noffdiag **2     !FZ: commented
      DO is = 1, nspin
        IF (ndiag .GT. 0) THEN
          DO ib = diag_nmin, diag_nmax
            WRITE (unit, 102) is, ib, mtxeld &
              (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)

  RETURN

  101 FORMAT (3F13.9, 1I8)   !FZ: changed from 2I8 to 1I8
  102 FORMAT (2I8, 2F15.9)
  103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_ekin_g

!-------------------------------------------------------------------------------

SUBROUTINE write_vnl (output_file_name, diag_nmin, diag_nmax, &     !FZ: all function is for metaGGA
  offdiag_nmin, offdiag_nmax)                     
!FZ: output all nonlocal energy for each bands and kpts (only diagonal)

  USE constants, ONLY : rytoev
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc, ehart !FZ:  added ehart
  USE exx, ONLY : vexx
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  !USE funct, ONLY : exx_is_active      !FZ: commented
  USE funct, ONLY : exx_is_active, dft_is_meta, get_meta     !FZ: Added for metaGGA commented in v_h
  USE gvect, ONLY : ngm, g
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, v_of_0, kedtau !, rho_core, rhog_core  FZ: for vhartree   !FZ: v_of_0 is the V(G=0) term, need to add to diagonal term which is g2kin,  because in v_H calculation, we did not calculate V(G=0)
  USE wavefunctions, ONLY : evc, psic
  USE wvfct, ONLY : npwx, nbnd, g2kin   !FZ: added g2kin
  USE g_psi_mod,            ONLY : h_diag, s_diag !FZ: for nonlocal pseudo
  USE noncollin_module, ONLY: noncolin, npol  !FZ: for nonlocal pseudo
  USE uspp,                 ONLY : vkb, nkb, okvan, deeq, qq_at, qq_so, deeq_nc, indv_ijkb0 !FZ: test for nonlocal
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp !FZ: test for nonlocal
  USE wvfct, ONLY: npwx, current_k !FZ: test for nonlocal, current_k is added for the use of h_psi_meta
  USE lsda_mod, ONLY: current_spin !FZ: test for nonlocal
  USE uspp_param, ONLY: upf, nh !FZ: test for nonlocal
  USE spin_orb, ONLY: lspinorb !FZ: test for nonlocal
  USE becmod,   ONLY : bec_type, becp, calbec, &  !FZ: added for nonlocal energies
                         allocate_bec_type, deallocate_bec_type  !FZ: added for nonlocal energies

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax
  integer, intent (inout) :: offdiag_nmin
  integer, intent (inout) :: offdiag_nmax

  integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, ib2, ikk
  integer, external :: global_kpoint_index
  complex (DP) :: dummy
  complex (DP), allocatable :: mtxeld (:, :)
  real (DP), allocatable :: vxcr (:, :)
  !real (DP), allocatable :: v_h (:, :)       !FZ: for metaGGA , stores hartree potential
  !real (DP), allocatable :: kedtaur (:, :)   !FZ: for metaGGA   kedtau: local K energy density,  kedtaur (in realspace)
  complex (DP), allocatable :: psic2 (:)
  complex (DP), allocatable :: hpsi (:)
  !REAL(DP)  ::  charge       !FZ: for metaGGA, stores total hartree energy, charge
  character(len=20) :: Ctemp
  INTEGER :: ierr   !FZ: for nonlocal
  !INTEGER, INTENT(IN)  :: lda, n, m
  !COMPLEX(DP), INTENT(INOUT) :: hpsi(lda*npol,m)
  !FZ: here we need lda, npol, m variables , but note that npwx is lda, we can use npwx as lda instead
  !FZ: npol is included from "USE noncollin_module, ONLY: noncolin, npol" 
  !FZ: n is number of plane waves in specific k point, it is npw, we specify it later npw = ngk ( ik - iks + 1 )
  !FZ: m is number of band, here we can use nbnd , (diag_nmax is nband) 
  !FZ: nbnd is the total number of bands used in nscf calculation (for Si it is 33)
  !COMPLEX(DP), INTENT(OUT) :: hpsinl(npwx*npol,nbnd)   !FZ: the hpsi to collect nonlocal potential (and a special potential in metaGGA (h_psi_meta)) 
  !COMPLEX(DP), INTENT(OUT) :: psinl(npwx*npol,nbnd)   !FZ: the psi used to calculate becp (we need  CALL calbec ( n, vkb, psi_nl, becp, m )) 
  COMPLEX(DP), allocatable :: hpsinl(:,:)   !FZ: the hpsi to collect nonlocal potential (and a special potential in metaGGA (h_psi_meta)) 
  COMPLEX(DP), allocatable :: psinl(:,:)   !FZ: the psi used to calculate becp (we need  CALL calbec ( n, vkb, psi_nl, becp, m )) 

  if(diag_nmin > diag_nmax) then
    call errore ( 'write_v_h_g', 'diag_nmin > diag_nmax', diag_nmin )  !FZ:
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    diag_nmax = nbnd
  ENDIF
  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)
  
  IF (ndiag .EQ. 0) RETURN

  unit = 4

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = (0.0D0, 0.0D0)
  ENDIF

  ALLOCATE (vxcr (dfftp%nnr, nspin))
  !ALLOCATE (v_h (dfftp%nnr, nspin))     !FZ: for meta GGA
  ALLOCATE (hpsi (dfftp%nnr))
  ALLOCATE (hpsinl (npwx*npol, nbnd))
  ALLOCATE (psinl (npwx*npol, nbnd))
  !ALLOCATE (kedtaur (dfftp%nnr, nspin))   !FZ: for meta GGA
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm )   !FZ: important! initialize becp!

  vxcr (:, :) = 0.0D0
  hpsinl(:,:) = (0.0D0, 0.0D0)  !FZ: for nonlocal potential , initialize
  psinl(:,:) = (0.0D0, 0.0D0)  !FZ: for nonlocal potential , initialize
  !v_h (:, :) = 0.0D0        !FZ: for metaGGA
  !kedtaur (:, :) = 0.0D0    !FZ: for metaGGA
  ehart  = 0.D0
  !charge = 0.D0
  !CALL v_h_only (rho%of_g(:,1), ehart, charge, v_h)   !FZ: for metaGGA   important check!rho%of_g(:,1)
  ALLOCATE( h_diag( npwx, npol ), STAT=ierr )    !FZ: for nonlocal
  IF( ierr /= 0 ) &    !FZ: for nonlocal
     CALL errore( ' diag_bands ', ' cannot allocate h_diag ', ABS(ierr) )     !FZ: for nonlocal
  ALLOCATE( s_diag( npwx, npol ), STAT=ierr )     !FZ: for nonlocal
  IF( ierr /= 0 ) &     !FZ: for nonlocal
     CALL errore( ' diag_bands ', ' cannot allocate s_diag ', ABS(ierr) )     !FZ: for nonlocal
  h_diag (:,:) = 0.0D0  !FZ: for nonlocal
  s_diag (:,:) = 0.0D0  !FZ: for nonlocal
  IF ( ionode ) WRITE (6, *) "nkb = ", nkb  !FZ: for future debug, when nkb > 0, we should calculate CALL add_vuspsi( lda, n, m, hpsi )
  
  !IF (ionode) THEN                                          !FZ: test
  !  OPEN (unit = unit, file = 'test_enl_h_diag_output', &          !FZ: test
  !    form = 'formatted', status = 'replace')               !FZ: test
  !    WRITE (unit, *) "h_diag = ", h_diag  !FZ:  test
  !    WRITE (unit, *) "s_diag = ", s_diag  !FZ:  test
  !    WRITE (unit, *) "vkb = ", vkb  !FZ:  test
  !    WRITE (unit, *) "npol = ", npol  !FZ:  test
  !  CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  !ENDIF                                          !FZ:  test

  DO ik = iks, ike
    ikk = ik - iks + 1
    npw = ngk ( ik - iks + 1 )
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    CALL threaded_memcpy(psinl, evc, nbnd*npol*npwx*2)
    vkb (:,:) = 0.0D0  !FZ: for nonlocal
    IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb ) !FZ: for nonlocal   !important : initialize vkb!
    h_diag (:,:) = 0.0D0  !FZ: for nonlocal
    s_diag (:,:) = 0.0D0  !FZ: for nonlocal
    hpsinl (:,:) = 0.0D0  !FZ: for nonlocal
    CALL usnldiag( npw, h_diag, s_diag )   !FZ: add nonlocal pseudopotential term to diagonal part of Hamiltonian
    g2kin(:) = 0.0D0  !FZ: important: initialize g2kin
    call g2_kin( ik )    !FZ:  check: if index is ik or ikk
    IF (abs(g2kin(1)) .lt. 0.0001_dp) THEN     !FZ: test
    IF (ionode) THEN                                          !FZ: test
      OPEN (unit = unit, file = 'test_enl_psi_nl_output', &          !FZ: test
        form = 'formatted', status = 'replace')               !FZ: test
        WRITE (unit, *) "psi_nl = ", psinl  !FZ:  test
        WRITE (unit, *) "npol = ", npol  !FZ:  test
        WRITE (unit, *) "nbnd = ", nbnd  !FZ:  test
        !WRITE (unit, *) "becp = ", becp  !FZ:  test
      CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
    ENDIF                                          !FZ:  test
    ENDIF                                          !FZ:  test
    CALL calbec ( npw, vkb, psinl, becp, nbnd )   !FZ: becp  ! <beta|psi>    original form: calbec ( n, vkb, psi, becp, m ) 
    !FZ: Here vkb is input, becp is output, so we need to initialize vkb before . vkb depends on k, so for each k point, we need to call initus_2
    CALL add_vuspsi( npwx, npw, nbnd, hpsinl )   !FZ:  CALL add_vuspsi( lda, n, m, hpsi )
    current_k = ik      !FZ: for metaGGA important, because h_psi_meta uses the variable current_k, we need to specify it here
    if (dft_is_meta()) call h_psi_meta (npwx, npw, nbnd, psinl, hpsinl)   !FZ: add specific contribution from metaGGA to hpsi (when call this, it is added to hpsi)   
    !FZ: original h_psi_meta (lda, n, m, psi, hpsi)
    IF (ionode) THEN              !FZ: test
      write(Ctemp,"(I2)") ik                                          !FZ: test
      OPEN (unit = unit, file = 'test_enl_output' // Trim(AdjustL(Ctemp)) //'.dat' , &  !FZ: test
        form = 'formatted', status = 'replace')               !FZ: test
        WRITE (unit, *) "hpsinl = ", hpsinl      !FZ test
        WRITE (unit, *) "vkb = ", vkb  !FZ:  test
      CLOSE (unit = unit, status = 'keep')      !FZ:  test
    ENDIF                                       !FZ:  test
    IF (abs(g2kin(1)) .lt. 0.0001_dp) THEN     !FZ: test
    IF (ionode) THEN                                          !FZ: test
      OPEN (unit = unit, file = 'test_enl_h_diag_output2', &          !FZ: test
        form = 'formatted', status = 'replace')               !FZ: test
        WRITE (unit, *) "h_diag = ", h_diag  !FZ:  test
        WRITE (unit, *) "s_diag = ", s_diag  !FZ:  test
      CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
    ENDIF                                          !FZ:  test
    ENDIF                                          !FZ:  test
    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        psic (:) = (0.0D0, 0.0D0)
        hpsi (:) = (0.0D0, 0.0D0)  !FZ: added initialization
        DO ig = 1, npw
          psic (igk_k(ig,ikk)) = evc (ig, ib)      !FZ:   evc stores all eigenfunctions in G space, use davcio function to extract ik - iks + 1 th kpoint 's all the eigen functions for all bands, 
        ENDDO
        !CALL invfft ('Rho', psic, dfftp)
        !DO ir = 1, dfftp%nnr
        !  psic (ir) = psic (ir) * v_h (ir, isk (ik))      !FZ:   v_hartree | psi>
        !ENDDO
        !CALL fwfft ('Rho', psic, dfftp)
        !hpsi (:) = (0.0D0, 0.0D0)
        !DO ig = 1, npw
        !  hpsi (ig) = psic (dfftp%nl (igk_k(ig,ikk)))
        !ENDDO
        DO ig = 1, npw                                  !FZ: calculate |k+G|^2 |psi>
          hpsi (ig) = g2kin(ig) * psic (igk_k(ig,ikk))   !FZ: calculate |k+G|^2 |psi>
          !hpsi (ig) = (g2kin(ig) + h_diag(ig, 1)) * psic (igk_k(ig,ikk))   !FZ: calculate |k+G|^2 |psi>
        ENDDO                                           !FZ: calculate |k+G|^2 |psi>
        !hpsi(1) = hpsi(1) + v_of_0 * psic (igk_k(1,ikk))  !FZ: add local potential 
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (ig) = evc (ig, ib)      !FZ: 
          !psic (igk_k(ig,ikk)) = evc (ig, ib)  !FZ: wrong
        ENDDO
        !IF (exx_is_active ()) &
        !   CALL vexx (npwx, npw, 1, psic, hpsi)
        dummy = (0.0D0, 0.0D0)
        DO ig = 1, npw
          !dummy = dummy + conjg (psic (ig)) * hpsi (ig)
          dummy = dummy + conjg (psic (ig)) * hpsinl (ig, ib)
        ENDDO
        dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
        CALL mp_sum (dummy, intra_bgrp_comm)
        mtxeld (ib - diag_nmin + 1, ik) = dummy
      ENDDO
    ENDIF
    
  ENDDO

  DEALLOCATE (vxcr)
  !DEALLOCATE (v_h)          !FZ: for metaGGA calculate hartree
  !DEALLOCATE (kedtaur)      !FZ: for metaGGA
  DEALLOCATE (hpsi)
  DEALLOCATE (hpsinl) !FZ: for nonlocal potential
  DEALLOCATE (psinl)  !FZ: for nonlocal potential
  DEALLOCATE (h_diag)  !FZ: for nonlocal
  DEALLOCATE (s_diag)  !FZ: for nonlocal
  CALL deallocate_bec_type ( becp )   !FZ: important! deallocate becp

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    DO ik = 1, nkstot / nspin
      WRITE (unit, 101) xk(:, ik), nspin * ndiag!, &   !FZ: commented
        !nspin * noffdiag **2     !FZ: commented
      DO is = 1, nspin
        IF (ndiag .GT. 0) THEN
          DO ib = diag_nmin, diag_nmax
            WRITE (unit, 102) is, ib, mtxeld &
              (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)

  RETURN

  101 FORMAT (3F13.9, 1I8)   !FZ: changed from 2I8 to 1I8
  102 FORMAT (2I8, 2F15.9)
  103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_vnl

!-------------------------------------------------------------------------------

SUBROUTINE write_vxc_g (output_file_name, diag_nmin, diag_nmax, &
  offdiag_nmin, offdiag_nmax, vxc_zero_rho_core)

  USE constants, ONLY : rytoev
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc
  USE exx, ONLY : vexx
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  !USE funct, ONLY : exx_is_active      !FZ: commented
  USE funct, ONLY : exx_is_active, dft_is_meta, get_meta     !FZ: Added for metaGGA
  USE gvect, ONLY : ngm, g
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE wavefunctions, ONLY : evc, psic
  USE wvfct, ONLY : npwx, nbnd

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax
  integer, intent (inout) :: offdiag_nmin      
  integer, intent (inout) :: offdiag_nmax      
  logical, intent (in) :: vxc_zero_rho_core

  integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2, ikk
  integer, external :: global_kpoint_index
  complex (DP) :: dummy
  complex (DP), allocatable :: mtxeld (:, :)
  complex (DP), allocatable :: mtxelo (:, :, :)
  real (DP), allocatable :: vxcr (:, :)
  real (DP), allocatable :: kedtaur (:, :)   !FZ: for metaGGA   kedtau: local K energy density,  kedtaur (in realspace)
  complex (DP), allocatable :: psic2 (:)
  complex (DP), allocatable :: hpsi (:)

  if(diag_nmin > diag_nmax) then
    call errore ( 'write_vxc_g', 'diag_nmin > diag_nmax', diag_nmin )
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    diag_nmax = nbnd
  ENDIF
  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

  if(offdiag_nmin > offdiag_nmax) then
    call errore ( 'write_vxc_g', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
  endif
  IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
  IF (offdiag_nmax .GT. nbnd)  then
    write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
    offdiag_nmax = nbnd
  ENDIF
  noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

  IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

  unit = 4

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = (0.0D0, 0.0D0)
  ENDIF
  IF (noffdiag .GT. 0) THEN
    ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
    mtxelo (:, :, :) = (0.0D0, 0.0D0)
  ENDIF

  ALLOCATE (vxcr (dfftp%nnr, nspin))
  IF (noffdiag .GT. 0) ALLOCATE (psic2 (dfftp%nnr))
  ALLOCATE (hpsi (dfftp%nnr))

  ALLOCATE (kedtaur (dfftp%nnr, nspin))   !FZ: for meta GGA
  vxcr (:, :) = 0.0D0
  kedtaur (:, :) = 0.0D0    !FZ: for metaGGA
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  !
  IF ( ionode ) WRITE (stdout, *) "dft_is_meta() = ", dft_is_meta(), "get_meta() = ", get_meta()    !FZ: test
  IF (dft_is_meta() .and. (get_meta() /= 4)) then         !FZ: for metaGGA
     CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, vxcr, kedtaur )    !FZ: for metaGGA
  ELSE                            !FZ: for metaGGA
     CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxcr)
  ENDIF                           !FZ: for metaGGA
  !
  IF ( ionode ) WRITE (stdout, '(" Test that stdout write_vxc_g is called")')    !FZ: test
  IF ( ionode ) WRITE (6, '(" Test that WRITE(6, ...) write_vxc_g is called")')  !FZ: test
  !IF (ionode) THEN                                          !FZ: test
  !  OPEN (unit = unit, file = 'test_pw2bgw_output_file', &  !FZ: test
  !    form = 'formatted', status = 'replace')               !FZ: test
  !    WRITE (unit, '(" Test output can be printed")') !FZ: test
  !    WRITE (unit, *) iks, ike    !FZ: test 
  IF (ionode) THEN                                          !FZ: test
    OPEN (unit = unit, file = 'test_vxcr_output', &  !FZ: test
      form = 'formatted', status = 'replace')               !FZ: test
      WRITE (unit, *) "vxcr = ", vxcr  !FZ:  test
    CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  ENDIF                                       !FZ:  test

  DO ik = iks, ike
    ikk = ik - iks + 1
    npw = ngk ( ik - iks + 1 )
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (dfftp%nl (igk_k(ig,ikk))) = evc (ig, ib)      !FZ:   evc stores all eigenfunctions in G space, use davcio function to extract ik - iks + 1 th kpoint 's all the eigen functions for all bands, 
   !FZ:   e.g.  silicon 10 kpts, 22 plane waves, so ig can be 1-22, igk_k(ig, ikk) reorders # of plane waves for each k point ikk
   !FZ:  igk_k is a 22*10 matrix
        ENDDO
        CALL invfft ('Rho', psic, dfftp)
        IF (ionode) THEN                                          !FZ: test
          OPEN (unit = unit, file = 'test_pw2bgw_output_3', &  !FZ: test
            form = 'formatted', status = 'replace')               !FZ: test
            WRITE (unit, *) "igk_k = ",igk_k  !FZ:  test
            WRITE (unit, *) "dfftp%nl = ",dfftp%nl  !FZ:  test
            WRITE (unit, *) "psic (dfftp%nl (igk_k(ig,ikk))) = ",psic (dfftp%nl (igk_k(1,1)))  !FZ:  test
            WRITE (unit, *) "npw = ",npw  !FZ:  test
          CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
        ENDIF                                       !FZ:  test
        DO ir = 1, dfftp%nnr
          psic (ir) = psic (ir) * vxcr (ir, isk (ik))
        ENDDO
        CALL fwfft ('Rho', psic, dfftp)
        hpsi (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          hpsi (ig) = psic (dfftp%nl (igk_k(ig,ikk)))
        ENDDO
        IF (ionode) THEN                                          !FZ: test
          OPEN (unit = unit, file = 'test_pw2bgw_output_2', &  !FZ: test
            form = 'formatted', status = 'replace')               !FZ: test
            WRITE (unit, *) "dfftp%nnr = ", dfftp%nnr," isk (ik) = ", isk    !FZ test     isk has dimesion = Max number of k-points (npk)
            WRITE (unit, *) "vxcr (ir=1, isk (ik=1)) ", vxcr (1, 1) !FZ:  test
            WRITE (unit, *) "iks = ",iks," ike = ", ike, " ikk = ", ikk      !FZ test
            WRITE (unit, *) "psic (dfftp%nl (igk_k(ig,ikk))) = ",psic (dfftp%nl (igk_k(1,1)))  !FZ:  test
            WRITE (unit, *) "dfftp%nl (igk_k(ig,ikk)) = ", dfftp%nl (igk_k(1,1)) !FZ:  test
          CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
        ENDIF                                       !FZ:  test
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (ig) = evc (ig, ib)
        ENDDO
        IF (exx_is_active ()) &
           CALL vexx (npwx, npw, 1, psic, hpsi)
        dummy = (0.0D0, 0.0D0)
        DO ig = 1, npw
          dummy = dummy + conjg (psic (ig)) * hpsi (ig)
        ENDDO
        IF (ionode) THEN                                          !FZ: test
          OPEN (unit = unit, file = 'test_pw2bgw_output_file_intermediate', &  !FZ: test
            form = 'formatted', status = 'replace')               !FZ: test
            WRITE (unit, '(" Test output can be printed")') !FZ: test
            WRITE (unit, *) "iks = ",iks," ike = ", ike, " ikk = ", ikk      !FZ test
            WRITE (unit, *) "psic (dfftp%nl (igk_k(ig,ikk))) = ",psic (dfftp%nl (igk_k(1,1)))  !FZ:  test
            WRITE (unit, *) "dfftp%nl (igk_k(ig,ikk)) = ", dfftp%nl (igk_k(1,1)) !FZ:  test
            WRITE (unit, *) "exx_is_active = ",exx_is_active ()  !FZ:  test
            WRITE (unit, *) "npw = ",npw," psic(1) = ", psic(1), " hpsi(1) = ", hpsi(1)      !FZ test
            WRITE (unit, *) "is = ",is," ib = ", ib, " dummy = ",dummy      !FZ test
          CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
        ENDIF                                       !FZ:  test
        dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
        CALL mp_sum (dummy, intra_bgrp_comm)
        mtxeld (ib - diag_nmin + 1, ik) = dummy
      ENDDO
      !WRITE (unit, '(" Test output can be printed")') !FZ: test
      !WRITE (unit, *) is, ib, dummy      !FZ test
      !WRITE (unit, 102) is, ib, dummy, dummy    !FZ  test
    ENDIF
    IF (noffdiag .GT. 0) THEN
      DO ib = offdiag_nmin, offdiag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (dfftp%nl (igk_k(ig,ikk))) = evc (ig, ib)
        ENDDO
        CALL invfft ('Rho', psic, dfftp)
        DO ir = 1, dfftp%nnr
          psic (ir) = psic (ir) * vxcr (ir, isk (ik))
        ENDDO
        CALL fwfft ('Rho', psic, dfftp)
        hpsi (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          hpsi (ig) = psic (dfftp%nl (igk_k (ig,ikk)))
        ENDDO
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (ig) = evc (ig, ib)
        ENDDO
        IF (exx_is_active ()) &
           CALL vexx (npwx, npw, 1, psic, hpsi)
        DO ib2 = offdiag_nmin, offdiag_nmax
          psic2 (:) = (0.0D0, 0.0D0)
          DO ig = 1, npw
            psic2 (ig) = evc (ig, ib2)
          ENDDO
          dummy = (0.0D0, 0.0D0)
          DO ig = 1, npw
            dummy = dummy + conjg (psic2 (ig)) * hpsi (ig)
          ENDDO
          dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
          CALL mp_sum (dummy, intra_bgrp_comm)
          mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin &
            + 1, ik) = dummy
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  !  CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  !ENDIF                                       !FZ:  test
  IF (ionode) THEN                                          !FZ: test
    OPEN (unit = unit, file = 'test_pw2bgw_output_mtxeld', &  !FZ: test
      form = 'formatted', status = 'replace')               !FZ: test
      WRITE (unit, *) "mtxeld = ", mtxeld  !FZ:  test
      WRITE (unit, *) "exx_is_active = ",exx_is_active ()  !FZ:  test
    CLOSE (unit = unit) !, status = 'keep')      !FZ:  test
  ENDIF                                       !FZ:  test

  DEALLOCATE (vxcr)
  DEALLOCATE (kedtaur)      !FZ: for metaGGA
  IF (noffdiag .GT. 0) DEALLOCATE (psic2)
  DEALLOCATE (hpsi)

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
  IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    DO ik = 1, nkstot / nspin
      WRITE (unit, 101) xk(:, ik), nspin * ndiag, &
        nspin * noffdiag **2
      DO is = 1, nspin
        IF (ndiag .GT. 0) THEN
          DO ib = diag_nmin, diag_nmax
            WRITE (unit, 102) is, ib, mtxeld &
              (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin)
          ENDDO
        ENDIF
        IF (noffdiag .GT. 0) THEN
          DO ib = offdiag_nmin, offdiag_nmax
            DO ib2 = offdiag_nmin, offdiag_nmax
              WRITE (unit, 103) is, ib2, ib, mtxelo &
                (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                ik + (is - 1) * nkstot / nspin)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
  IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

  RETURN

  101 FORMAT (3F13.9, 2I8)
  102 FORMAT (2I8, 2F15.9)
  103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_vxc_g

!-------------------------------------------------------------------------------

SUBROUTINE write_vscg ( output_file_name, real_or_complex, symm_type )

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  USE constants, ONLY : pi, tpi, eps6
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, ecutrho
  USE io_global, ONLY : ionode
  USE ions_base, ONLY : nat, atm, ityp, tau 
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : vltot, v
  USE symm_base, ONLY : s, ftau, nsym
  USE wavefunctions, ONLY : psic
  USE matrix_inversion

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  integer, intent (in) :: real_or_complex
  character ( len = 9 ), intent (in) :: symm_type

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  integer :: unit, id, is, ir, ig, i, j, k, ierr
  integer :: nd, ns, nr, ng_l, ng_g
  integer :: ntran, cell_symmetry, nrecord
  real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: g_g ( :, : )
  real (DP), allocatable :: vscr_g ( :, : )
  complex (DP), allocatable :: vscg_g ( :, : )

  INTEGER, EXTERNAL :: atomic_number

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  ! this is supposed to be VSC-Real/Complex but BGW wfn_rho_vxc IO
  ! does not recognize VSC header so we are using VXC instead
  IF ( real_or_complex .EQ. 1 ) THEN
    WRITE ( stitle, '("VXC-Real",24X)' )
  ELSE
    WRITE ( stitle, '("VXC-Complex",21X)' )
  ENDIF

  unit = 4
  nrecord = 1
  nd = 3

  ns = nspin
  nr = dfftp%nnr
  ng_l = ngm
  ng_g = ngm_g

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_vscg', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2 )
    t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( dfftp%nr1 )
    t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( dfftp%nr2 )
    t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( dfftp%nr3 )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  ALLOCATE ( g_g ( nd, ng_g ) )
  ALLOCATE ( vscr_g ( ng_g, ns ) )
  ALLOCATE ( vscg_g ( ng_g, ns ) )

  DO ig = 1, ng_g
    DO id = 1, nd
      g_g ( id, ig ) = 0
    ENDDO
  ENDDO
  DO is = 1, ns
    DO ig = 1, ng_g
      vscg_g ( ig, is ) = ( 0.0D0, 0.0D0 )
    ENDDO
  ENDDO

  DO ig = 1, ng_l
    g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  vscr_g ( :, : ) = 0.0D0
  DO is = 1, ns
    DO ir = 1, nr
      psic ( ir ) = CMPLX ( v%of_r ( ir, is ) + vltot ( ir ), 0.0D0, KIND=dp )   !FZ: vsc is v%of_r + vltot
    ENDDO
    CALL fwfft ( 'Rho', psic, dfftp )
    DO ig = 1, ng_l
      vscg_g ( ig_l2g ( ig ), is ) = psic ( dfftp%nl ( ig ) )
    ENDDO
  ENDDO

  CALL mp_sum ( g_g, intra_bgrp_comm )
  CALL mp_sum ( vscg_g, intra_bgrp_comm )

  IF ( ionode ) THEN
    OPEN ( unit = unit, file = TRIM ( output_file_name ), &
      form = 'unformatted', status = 'replace' )
    WRITE ( unit ) stitle, sdate, stime
    WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho
    WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    IF ( real_or_complex .EQ. 1 ) THEN
      WRITE ( unit ) ( ( dble ( vscg_g ( ig, is ) ), &
        ig = 1, ng_g ), is = 1, ns )
    ELSE
      WRITE ( unit ) ( ( vscg_g ( ig, is ), &
        ig = 1, ng_g ), is = 1, ns )
    ENDIF
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( vscg_g )
  DEALLOCATE ( vscr_g )
  DEALLOCATE ( g_g )

  RETURN

END SUBROUTINE write_vscg

!-------------------------------------------------------------------------------

SUBROUTINE write_vkbg (output_file_name, symm_type, wfng_kgrid, &
  wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3)

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  USE constants, ONLY : pi, tpi, eps6
  USE fft_base, ONLY : dfftp
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, g, mill, ecutrho
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat, atm, ityp, tau, nsp
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, wk, ngk, nks, nkstot, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum, mp_max, mp_get, mp_barrier
  USE mp_world, ONLY : mpime, nproc, world_comm
  USE mp_bands, ONLY : intra_bgrp_comm, nbgrp
  USE mp_pools, ONLY : me_pool, root_pool, npool, nproc_pool, &
    intra_pool_comm, inter_pool_comm
  USE mp_wave, ONLY : mergewf
  USE start_k, ONLY : nk1, nk2, nk3, k1, k2, k3
  USE symm_base, ONLY : s, ftau, nsym
  USE uspp, ONLY : nkb, vkb, deeq
  USE uspp_param, ONLY : nhm, nh
  USE wvfct, ONLY : npwx, nbnd, g2kin   !FZ: added nbnd g2kin
  USE gvecw, ONLY : ecutwfc
  USE matrix_inversion

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  character ( len = 9 ), intent (in) :: symm_type
  logical, intent (in) :: wfng_kgrid
  integer, intent (in) :: wfng_nk1
  integer, intent (in) :: wfng_nk2
  integer, intent (in) :: wfng_nk3
  real (DP), intent (in) :: wfng_dk1
  real (DP), intent (in) :: wfng_dk2
  real (DP), intent (in) :: wfng_dk3

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  integer :: i, j, k, ierr, ik, is, ig, ikb, iat, isp, ih, jh, &
    unit, iks, ike, npw, npw_g, npwx_g, ngg, ipsour, &
    igwx, local_pw, id, nd, ntran, cell_symmetry, nrecord
  real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: kmap ( : )
  integer, allocatable :: smap ( : )
  integer, allocatable :: gvec ( :, : )
  integer, allocatable :: ngk_g ( : )
  integer, allocatable :: igk_l2g ( :, : )
  integer, allocatable :: itmp ( : )
  integer, allocatable :: igwk ( : )
  integer, allocatable :: igwf_l2g ( : )
  integer, allocatable :: ipmask ( : )
  complex (DP), allocatable :: vkb_g ( : )

  INTEGER, EXTERNAL :: atomic_number, global_kpoint_index

  IF ( nkb == 0 ) RETURN

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  ! BGW wfn_rho_vxc IO does not recognize VKB header so this file
  ! is read directly by SAPO code in BerkeleyGW
  WRITE ( stitle, '("VKB-Complex",21X)' )

  unit = 4
  nrecord = 1
  nd = 3

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_vkbg', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2 )
    t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( dfftp%nr1 )
    t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( dfftp%nr2 )
    t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( dfftp%nr3 )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  ALLOCATE ( kmap ( nkstot ) )
  ALLOCATE ( smap ( nkstot ) )

  DO i = 1, nkstot
    j = ( i - 1 ) / nspin
    k = i - 1 - j * nspin
    kmap ( i ) = j + k * ( nkstot / nspin ) + 1
    smap ( i ) = k + 1
  ENDDO
  ierr = 0
  DO i = 1, nkstot
    ik = kmap ( i )
    is = smap ( i )
    IF ( ik .GE. iks .AND. ik .LE. ike .AND. is .NE. isk ( ik ) ) &
      ierr = ierr + 1
  ENDDO
  CALL mp_max ( ierr, world_comm )
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_vkbg', 'smap', ierr )

  ALLOCATE ( gvec ( 3, ngm_g ) )
  gvec = 0
  DO ig = 1, ngm
    gvec ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    gvec ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    gvec ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  CALL mp_sum ( gvec, intra_bgrp_comm )

  ALLOCATE ( ngk_g ( nkstot ) )
  ALLOCATE ( igk_l2g ( npwx, nks ) )
  ngk_g = 0
  igk_l2g = 0
  DO ik = 1, nks
    npw = ngk ( ik )
    DO ig = 1, npw
      igk_l2g ( ig, ik ) = ig_l2g ( igk_k ( ig, ik ) )
    ENDDO
  ENDDO
  DO ik = 1, nks
    ngk_g ( ik + iks - 1 ) = ngk ( ik )
  ENDDO
  CALL mp_sum ( ngk_g, inter_pool_comm )
  CALL mp_sum ( ngk_g, intra_pool_comm )
  ngk_g = ngk_g / nbgrp
  npw_g = MAXVAL ( igk_l2g ( :, : ) )
  CALL mp_max ( npw_g, intra_pool_comm )
  npwx_g = MAXVAL ( ngk_g ( : ) )

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'unformatted', status = 'replace')
    WRITE ( unit ) stitle, sdate, stime
    WRITE ( unit ) nspin, ngm_g, ntran, cell_symmetry, nat, ecutrho, &
      nkstot / nspin, nsp, nkb, nhm, npwx_g, ecutwfc
    IF ( wfng_kgrid ) THEN
      WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, wfng_nk1, wfng_nk2, wfng_nk3, &
        wfng_dk1, wfng_dk2, wfng_dk3
    ELSE
      WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, nk1, nk2, nk3, &
        0.5D0 * dble ( k1 ), 0.5D0 * dble ( k2 ), 0.5D0 * dble ( k3 )
    ENDIF
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) ( ngk_g ( ik ), ik = 1, nkstot / nspin )
    WRITE ( unit ) ( wk ( ik ) * dble ( nspin ) / 2.0D0, ik = 1, nkstot / nspin )
    WRITE ( unit ) ( ( xk ( id, ik ), id = 1, nd ), ik = 1, nkstot / nspin )
    WRITE ( unit ) ( ityp ( iat ), iat = 1, nat )
    WRITE ( unit ) ( nh ( isp ), isp = 1, nsp )
    WRITE ( unit ) ( ( ( ( deeq ( jh, ih, iat, is ), &
      jh = 1, nhm ), ih = 1, nhm ), iat = 1, nat ), is = 1, nspin )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ngm_g
    WRITE ( unit ) ( ( gvec ( id, ig ), id = 1, nd ), ig = 1, ngm_g )
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  ALLOCATE ( igwk ( npwx_g ) )

  DO i = 1, nkstot

    ik = kmap ( i )
    is = smap ( i )

    igwk = 0

    ALLOCATE ( itmp ( npw_g ) )
    itmp = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
      DO ig = 1, ngk ( ik - iks + 1 )
        itmp ( igk_l2g ( ig, ik - iks + 1 ) ) = igk_l2g ( ig, ik - iks + 1 )
      ENDDO
    ENDIF
    CALL mp_sum ( itmp, intra_bgrp_comm )
    ngg = 0
    DO ig = 1, npw_g
      IF ( itmp ( ig ) .EQ. ig ) THEN
        ngg = ngg + 1
        igwk ( ngg ) = ig
      ENDIF
    ENDDO
    DEALLOCATE ( itmp )

    IF ( ionode ) THEN
      IF ( is .EQ. 1 ) THEN
        WRITE ( unit ) nrecord
        WRITE ( unit ) ngk_g ( ik )
        WRITE ( unit ) ( ( gvec ( j, igwk ( ig ) ), j = 1, 3 ), &
          ig = 1, ngk_g ( ik ) )
      ENDIF
    ENDIF

    local_pw = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
      npw = ngk ( ik - iks + 1 )
      CALL init_us_2 ( npw, igk_k(1, ik-iks+1), xk ( 1, ik ), vkb )
      local_pw = npw
    ENDIF

    ALLOCATE ( igwf_l2g ( local_pw ) )
    igwf_l2g = 0
    DO ig = 1, local_pw
      ngg = igk_l2g ( ig, ik - iks + 1 )
      DO j = 1, ngk_g ( ik )
        IF ( ngg .EQ. igwk ( j ) ) THEN
          igwf_l2g ( ig ) = j
          EXIT
        ENDIF
      ENDDO
    ENDDO

    ALLOCATE ( ipmask ( nproc ) )
    ipmask = 0
    ipsour = ionode_id
    IF ( npool .GT. 1 ) THEN
      IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
        IF ( me_pool .EQ. root_pool ) ipmask ( mpime + 1 ) = 1
      ENDIF
      CALL mp_sum ( ipmask, world_comm )
      DO j = 1, nproc
        IF ( ipmask ( j ) .EQ. 1 ) ipsour = j - 1
      ENDDO
    ENDIF
    DEALLOCATE ( ipmask )

    igwx = 0
    ierr = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) &
      igwx = MAXVAL ( igwf_l2g ( 1 : local_pw ) )
    CALL mp_max ( igwx, intra_pool_comm )
    IF ( ipsour .NE. ionode_id ) &
      CALL mp_get ( igwx, igwx, mpime, ionode_id, ipsour, 1, world_comm )
    ierr = 0
    IF ( ik .GE. iks .AND. ik .LE. ike .AND. igwx .NE. ngk_g ( ik ) ) &
      ierr = 1
    CALL mp_max ( ierr, world_comm )
    IF ( ierr .GT. 0 ) &
      CALL errore ( 'write_vkbg', 'igwx ngk_g', ierr )

    ALLOCATE ( vkb_g ( MAX ( 1, igwx ) ) )

    DO ikb = 1, nkb

      vkb_g = ( 0.0D0, 0.0D0 )
      IF ( npool .GT. 1 ) THEN
        IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
          CALL mergewf ( vkb ( :, ikb ), vkb_g, local_pw, igwf_l2g, &
            me_pool, nproc_pool, root_pool, intra_pool_comm )
        ENDIF
        IF ( ipsour .NE. ionode_id ) THEN
          CALL mp_get ( vkb_g, vkb_g, mpime, ionode_id, ipsour, ikb, &
            world_comm )
        ENDIF
      ELSE
        CALL mergewf ( vkb ( :, ikb ), vkb_g, local_pw, igwf_l2g, &
          mpime, nproc, ionode_id, world_comm )
      ENDIF

      IF ( ionode ) THEN
        WRITE ( unit ) nrecord
        WRITE ( unit ) ngk_g ( ik )
        WRITE ( unit ) ( vkb_g ( ig ), ig = 1, igwx )
      ENDIF
      call g2_kin( ik )    !FZ:  check: if index is ik or ikk
      IF (abs(g2kin(1)) .lt. 0.0001_dp) THEN     !FZ: test
      IF (ionode) THEN                                          !FZ: test
        OPEN (unit = 4, file = 'test_pw2bgw_vkb', &  !FZ: test
          form = 'formatted', status = 'replace')               !FZ: test
          WRITE (unit, *) "vkb_g = ", vkb_g  !FZ:  test
          WRITE (unit, *) "ngk_g = ", ngk_g  !FZ:  test
          WRITE (unit, *) "vkb = ", vkb  !FZ:  test
        CLOSE (unit = 4) !, status = 'keep')      !FZ:  test
      ENDIF                                       !FZ:  test
      ENDIF                                       !FZ:  test

    ENDDO

    DEALLOCATE ( vkb_g )
    DEALLOCATE ( igwf_l2g )

  ENDDO

  IF ( ionode ) THEN
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( igwk )
  DEALLOCATE ( igk_l2g )
  DEALLOCATE ( ngk_g )
  DEALLOCATE ( gvec )
  DEALLOCATE ( smap )
  DEALLOCATE ( kmap )

  RETURN

END SUBROUTINE write_vkbg

!-------------------------------------------------------------------------------

subroutine check_inversion(real_or_complex, ntran, mtrx, nspin, warn, real_need_inv, tnp)

! check_inversion    Originally By David A. Strubbe    Last Modified 11/18/2013
! Check whether our choice of real/complex version is appropriate given the
! presence or absence of inversion symmetry.

  USE constants, ONLY : eps6
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP

  implicit none

  integer, intent(in) :: real_or_complex
  integer, intent(in) :: ntran
  integer, intent(in) :: mtrx(3, 3, 48) !< symmetry operations matrices
  integer, intent(in) :: nspin
  logical, intent(in) :: warn !< set to false to suppress warnings, for converters
  logical, intent(in) :: real_need_inv !< use for generating routines to block real without inversion
     !! this is not always true so that it is possible to run real without using symmetries
  real(DP), intent(in) :: tnp(3, 48) !< fractional translations.
     !! optional only to avoid changing external interface for library.

  integer :: invflag, isym, ii, jj, itest
  logical :: origin_inv

  invflag = 0
  origin_inv = .false.
  do isym = 1, ntran
    itest = 0
    do ii = 1, 3
      do jj = 1, 3
        if(ii .eq. jj) then
          itest = itest + (mtrx(ii, jj, isym) + 1)**2
        else
          itest = itest + mtrx(ii, jj, isym)**2
        endif
      enddo
    enddo
    if(itest .eq. 0) then
      invflag = invflag + 1
      !if(present(tnp)) then
        if(sum(abs(tnp(1:3, isym))) < eps6) origin_inv = .true.
      !else
      !  origin_inv = .true.
      !endif
    endif
  enddo
  if(invflag > 0 .and. .not. origin_inv) then
    write(0, '(a)') "WARNING: Inversion symmetry is present only with a fractional translation."
    write(0, '(a)') "Apply the translation so inversion is about the origin, to be able to use the real version."
  endif
  if(invflag .gt. 1) write(0, '(a)') "WARNING: More than one inversion symmetry operation is present."

!  if(invflag > 0 .and. .not. present(tnp)) then
!    write(0, '(a)') "WARNING: check_inversion did not receive fractional translations."
!    write(0, '(a)') "Cannot confirm that inversion symmetry is about the origin for use of real version."
!  endif

  if(real_or_complex .eq. 2) then
    if(origin_inv .and. warn .and. nspin == 1) then
      if(ionode) &
        write(0, '(a)') "WARNING: Inversion symmetry about the origin is present. The real version would be faster."
    endif
  else
    if(.not. origin_inv) then
      if(real_need_inv) then
        call errore("check_inversion", "The real version cannot be used without inversion symmetry about the origin.", -1)
      endif
      if(ionode) then
        write(0, '(a)') "WARNING: Inversion symmetry about the origin is absent in symmetries used to reduce k-grid."
        write(0, '(a)') "Be sure inversion about the origin is still a spatial symmetry, or you must use complex version instead."
      endif
    endif
    if(nspin > 1) then
      call errore("check_inversion", "Real version may only be used for spin-unpolarized calculations.", nspin)
    endif
  endif

  return

end subroutine check_inversion

!-------------------------------------------------------------------------------

END PROGRAM pw2bgw

