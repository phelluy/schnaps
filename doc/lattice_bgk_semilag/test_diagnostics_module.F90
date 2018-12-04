module test_diagnostics_module
use, intrinsic :: iso_c_binding

  implicit none
  !include 'fftw3.f03'
  

contains

  subroutine init_sol_Wave_smooth( &
    x, &
    t, &
    w, &
    VMAX )
    real*8, intent(in) :: x
    real*8, intent(in) :: t
    complex*16, intent(out) :: w(:)
    real*8, intent(in), optional :: VMAX
    complex*16 ::rho, u
    real*8 :: VMAX_val

    if(present(VMAX))then
      VMAX_val = VMAX
    else
      VMAX_val = -1._8   
   endif

   rho = exp(- 30._8 * x * x)
   u = 0.0
   
   w(1) = VMAX_val*VMAX_val* 0.5 *rho
   w(2) = rho * (1_8 - VMAX_val)
   w(3) = VMAX_val*VMAX_val* 0.5 *rho

 end subroutine init_sol_Wave_smooth

 subroutine init_sol_Wave_nonsmooth( &
    x, &
    t, &
    w, &
    VMAX )
    real*8, intent(in) :: x
    real*8, intent(in) :: t
    complex*16, intent(out) :: w(:)
    real*8, intent(in), optional :: VMAX
    complex*16 ::rho, u
    real*8 :: VMAX_val

    if(present(VMAX)) then
       VMAX_val = VMAX
    else
       VMAX_val = -1._8   
    endif
    
    if(abs(x) <0.1) then
       rho = 1
    else
       rho =0
    endif
    u = 0.0
   
    w(1) = VMAX_val*VMAX_val* 0.5 *rho
    w(2) = rho * (1_8 - VMAX_val)
    w(3) = VMAX_val*VMAX_val* 0.5 *rho
    
 end subroutine init_sol_Wave_nonsmooth

 subroutine init_sol_Euler_nonsmooth( &
    x, &
    t, &
    w, &
    VMAX, &
    CSON)
    real*8, intent(in) :: x
    real*8, intent(in) :: t
    complex*16, intent(out) :: w(:)
    real*8, intent(in), optional :: VMAX
    real*8, intent(in), optional :: CSON
    
    real*8 :: xi
    real*8 :: RL,UL,UR,RR,ratio
    
    ratio=CSON/VMAX
    RL= 2._8
    UL = -0._8
    RR = 1._8
    UR = 0._8
    xi = x 
    if (xi < 0._8) then
      w(1) = RL * UL * ((UL/VMAX) - 1._8) / 2._8 &
           + ratio * ratio * RL / 2._8
      w(2) = RL * (1._8 - (UL * UL + CSON * CSON)/(VMAX*VMAX))
      w(3) = RL * UL * ((UL/VMAX) + 1._8) / 2._8 &
           + ratio * ratio * RL / 2._8
    else
      w(1) = RR * UR * ((UR/VMAX) - 1._8) / 2._8 &
           + ratio * ratio * RR / 2._8
      w(2) = RR * (1._8 - (UR * UR + CSON * CSON)/(VMAX*VMAX))
      w(3) = RR * UR * ((UR/VMAX) + 1._8) / 2._8 &
           + ratio * ratio * RR / 2._8
    endif    
    
  end subroutine init_sol_Euler_nonsmooth

  subroutine init_sol_Euler_smooth( &
    x, &
    t, &
    w, &
    VMAX, &
    CSON )
    real*8, intent(in) :: x
    real*8, intent(in) :: t
    complex*16, intent(out) :: w(:)
    real*8, intent(in), optional :: VMAX
    real*8, intent(in), optional :: CSON
    complex*16 ::rho, u
    real*8 :: CSON_val, ratio

    if(present(CSON))then
      CSON_val = CSON
    else
     CSON_val = 1._8   
  endif

   rho = 1.0+exp(- 30._8 * x * x)
   u = 0.0

   ratio=CSON_val/VMAX
   
   w(1) = 0.5 * rho * u * ((u/VMAX) - 1._8) + ratio*ratio* 0.5 *rho
   w(2) = rho * (1_8 - (u * u + CSON_val* CSON_val)/(VMAX*VMAX))
   w(3) = 0.5 * rho * u * ((u/VMAX) + 1._8) + ratio*ratio* 0.5 *rho

 end subroutine init_sol_Euler_smooth

  subroutine compute_exact_sol( &
    x, &
    t, &
    w, &
    VMAX, &
    RL, &
    UL, &
    RR, &
    UR, &
    CSON)
    real*8, intent(in) :: x
    real*8, intent(in) :: t
    complex*16, intent(out) :: w(:)
    real*8, intent(in), optional :: VMAX
    real*8, intent(in), optional :: RL
    real*8, intent(in), optional :: UL
    real*8, intent(in), optional :: RR
    real*8, intent(in), optional :: UR
    real*8, intent(in), optional :: CSON
    
    real*8 :: xi
    real*8 :: VMAX_val
    real*8 :: RL_val
    real*8 :: UL_val
    real*8 :: RR_val
    real*8 :: UR_val
    real*8 :: CSON_val
    
    if(present(VMAX))then
      VMAX_val = VMAX
    else
      VMAX_val = -1._8   
    endif
    
    if(present(RL))then
      RL_val = RL
    else
      RL_val = 2._8   
    endif
    if(present(UL))then
      UL_val = UL
    else
      UL_val = -0._8   
    endif

    if(present(RR))then
      RR_val = RR
    else
      RR_val = 1._8   
    endif
    if(present(UR))then
      UR_val = UR
    else
      UR_val = -0._8   
    endif

    if(present(RL))then
      CSON_val = CSON
    else
      CSON_val = 0.6_8   
    endif
    
    
    
    !Riemann problem test case
    xi = x + VMAX_val * t
    if (xi < 0._8) then
      w(1) = RL_val * UL_val * (UL_val - 1._8) / 2._8 &
        + CSON_val * CSON_val * RL_val / 2._8
    else
      w(1) = RR_val * UR_val * (UR_val - 1._8) / 2._8 &
        + CSON_val * CSON_val * RR_val / 2._8
    endif    
    xi = x
    if (xi < 0._8) then
      w(2) = RL_val * (1._8 - UL_val * UL_val - CSON_val * CSON_val)
    else
      w(2) = RR_val * (1._8 - UR_val * UR_val - CSON_val * CSON_val)
    endif
    xi = x - VMAX_val * t
    if (xi < 0._8) then
      w(3) = RL_val * UL_val * (UL_val + 1._8) / 2._8 &
        + CSON_val * CSON_val * RL_val / 2._8
    else
      w(3) = RR_val * UR_val * (UR_val + 1._8) / 2._8 &
        + CSON_val * CSON_val * RR_val / 2._8
    endif
    
  end subroutine compute_exact_sol



  subroutine compute_solex_rhou_smooth( &
    x, &
    t, &
    VMAX, &
    rho, &
    u )
    real*8, intent(in) :: x
    real*8, intent(in) :: t
    complex*16, intent(out) :: rho,u
    real*8, intent(in), optional :: VMAX
    
    real*8 :: xi, w1,w2
    real*8 :: VMAX_val

    if(present(VMAX))then
      VMAX_val = VMAX
    else
      VMAX_val = -1._8   
    endif

    xi = x - VMAX_val * t
    w1 = exp(- 30._8 * xi * xi)
    
    xi = x + VMAX_val * t
    w2 = exp(- 30._8 * xi * xi)

    rho = 0.5 * (w1 + w2)
    u = (0.5 * (w1 - w2)) / rho
    
  end subroutine compute_solex_rhou_smooth

  subroutine compute_solex_rhou_nonsmooth( &
    x, &
    t, &
    VMAX, &
    rho, &
    u )
    real*8, intent(in) :: x
    real*8, intent(in) :: t
    complex*16, intent(out) :: rho,u
    real*8, intent(in), optional :: VMAX
    
    real*8 :: xi, w1,w2
    real*8 :: VMAX_val

    if(present(VMAX))then
      VMAX_val = VMAX
    else
      VMAX_val = -1._8   
    endif

    xi = x - VMAX_val * t
    if(abs(x) <0.1) then
       w1 = 1
    else
       w1 =0
    endif
    
    xi = x + VMAX_val * t
    if(abs(x) <0.1) then
       w2 = 1
    else
       w2 =0
    endif

    rho = 0.5 * (w1 + w2)
    u = (0.5 * (w1 - w2)) / rho
    
  end subroutine compute_solex_rhou_nonsmooth

   subroutine write_wave_solex( &
    file_id, &
    xmin, &
    xmax, &
    w, &
    VMAX, &
    smooth, &
    tmax )
    integer, intent(in) :: file_id
    real*8, intent(in) :: xmin
    real*8, intent(in) :: xmax
    real*8, intent(in) :: tmax
    real*8, intent(in) :: vmax
    complex*16, intent(in) :: w(:,:)
    logical, intent(in) :: smooth
    
    integer :: N
    integer :: M
    integer :: i
    real*8 :: dx
    complex*16 :: rho, rho_ex
    complex*16 :: u, u_ex
    real*8 :: x
    
    N = size(w,2)-1
    M = size(w,1)
    
    dx = (xmax-xmin)/real(N,8)
    
    do i=1,N+1
       x = xmin+real(i-1,8)*dx
       if(abs(x) .le. 1_8) then
          if(smooth) then
             call compute_solex_rhou_smooth(x,tmax,VMAX,rho_ex,u_ex)
          else
             call compute_solex_rhou_nonsmooth(x,tmax,VMAX,rho_ex,u_ex)
          endif
          write(file_id,'(5g25.12)') x,real(rho_ex,8),real(u_ex,8)
       endif
    enddo
    
  end subroutine write_wave_solex

  subroutine write_sol( &
    file_id, &
    xmin, &
    xmax, &
    VMAX, &
    w)
    integer, intent(in) :: file_id
    real*8, intent(in) :: xmin
    real*8, intent(in) :: xmax
    complex*16, intent(in) :: w(:,:)
    real*8, intent(in):: VMAX
    integer :: N
    integer :: M
    integer :: i,c
    real*8 :: dx
    complex*16 :: rho
    complex*16 :: u
    real*8 :: x
    
    N = size(w,2)-1
    M = size(w,1)
    
    dx = (xmax-xmin)/real(N,8)
    
    do i=1,N+1
       x = xmin+real(i-1,8)*dx 
       if((x .le. 2_8) .and. (x .ge. -2_8)) then
          c=c+1
          rho = sum(w(:,i))
          u = (VMAX*(w(M,i)-w(1,i)))/rho
          write(file_id,'(3g25.12)') x,real(rho,8),real(u,8)
       endif
    enddo
 
  end subroutine write_sol
 



  
end module test_diagnostics_module
