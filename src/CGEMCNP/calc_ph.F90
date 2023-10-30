!======================================================================     
  Subroutine calc_ph( ff,ff_new,dT, S, T, PAR, Wind, lat, d_sfc, is_bottom, is_day )

!======================================================================
  use cgem
  use mvars
  use schism_glbl, only : rkind

  IMPLICIT NONE

!---------------------------------------------
! Interface variables
!---------------------------------------------------------------------
  logical, intent(in)  :: is_bottom ! Is it the bottom?  For instant remineralization
  logical, intent(in)  :: is_day    ! Uptake only occurs during the day
  real(rkind), intent(in)     :: lat       ! For mocsy- latitude 
  real(rkind), intent(in)     :: d_sfc     ! For mocsy
  real(rkind), intent(in)     :: S,T,dT
  real(rkind),intent(in), dimension(nf) :: ff
  real(rkind),intent(out), dimension(nf) :: ff_new
!---------------------------------------------------------------------------------------
! Local Variables
!-----------------------------------------------------
  integer :: isp, isz ! Loop indicies, isp/isz is for phytoplankton/zooplankton species
!------------------------------------ 
! Phytoplankton parameters
! Phytoplankton uptake and growth
!Output vars for alkalinity subroutine:
  real(rkind),dimension(km) :: ph_calc, pco2_calc, fco2, co2, hco3, co3, omegaa, omegac, betad_calc 
  real(rkind) :: rhosw, p, tempis(1)
  real(rkind) :: patm = 1.
  real(rkind) :: m_alk, m_dic, m_si, m_po4
  real(rkind) :: m_lat,m_T,m_S,m_d_sfc
!mocsy needs lat to be an array
  m_lat = lat
  m_T = T
  m_S = S
  m_d_sfc = d_sfc

  !------------------------------------------------------------
  ! Carbon Chemistry
  !--------------------------------------------------------------
  !!! MOCSY alkalinity expressions:
  m_alk = ALK/1000.d0
  m_dic = DIC/1000.d0
  m_si  = Si/1000.d0
  m_po4 = PO4/1000.d0
  call vars(ph_calc, pco2_calc, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD_calc, rhoSW, p, tempis,&
  &    m_T, m_S, m_alk, m_dic, m_si, m_po4, patm, m_d_sfc, m_lat, 1, &
  &    'mol/m3', 'Tinsitu', 'm ', 'u74', 'l  ', 'pf ', 'Pzero  ')
  pH = ph_calc

  return
  end subroutine cgem_step
!---------------------------------------------------------------------- 
