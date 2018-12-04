module semi_lag_module
use, intrinsic :: iso_c_binding

  implicit none
!include 'fftw3.f03'

contains
  
  subroutine compute_dt( &
    vmax, &
    cfl, &
    dx, &
    c, &
    dt)
    
    real*8, intent(in) :: vmax
    real*8, intent(in) :: cfl
    real*8, intent(in) :: dx
    complex*16, intent(in) :: c
    complex*16, intent(out) :: dt    

    if(abs(vmax)>1e-10_8)then
      dt = cfl * dx / abs(vmax) * c
    else
      dt = cfl * dx * c
    endif


  end subroutine compute_dt
  

  subroutine compute_dt_big( &
    vmax, &
    cfl, &
    dx, &
    dt)
    
    real*8, intent(in) :: vmax
    real*8, intent(in) :: cfl
    real*8, intent(in) :: dx
    real*8, intent(out) :: dt    

    if(abs(vmax)>1e-10_8)then
      dt = cfl * dx / abs(vmax) 
    else
      dt = cfl * dx
    endif


  end subroutine compute_dt_big

  subroutine bgk_relax(w, baryc,VMAX,CSON,factor)
    complex*16, intent(inout) :: w(:)
    complex*16, intent(in) :: baryc
    real*8, intent(in) :: CSON
    real*8, intent(in) :: VMAX
    real*8, intent(in), optional :: factor
    complex*16 :: r
    complex*16 :: q
    complex*16 :: u
    complex*16 :: w_check(3)
    complex*16 :: w_eq(3)
    real*8 :: loc_factor, ratio
    
    if(present(factor))then
      loc_factor = factor
    else
      loc_factor = 1._8  
    endif
    
    r = w(1)+w(2)+w(3)
    q = VMAX*(w(3)-w(1))
    u = q/r

    ratio=CSON/VMAX
    
    w_eq(1) = (q /VMAX) * ( loc_factor * (u/VMAX) - 1._8) / 2._8 + ratio * ratio * r / 2._8
    w_eq(2) = r * (1._8 - (loc_factor*u * u + CSON * CSON)/(VMAX*VMAX))
    w_eq(3) = (q /VMAX) * ( loc_factor * (u/VMAX) + 1._8) / 2._8 + ratio * ratio* r / 2._8
    
    w(1) = w(1)*baryc+(1._8-baryc) * w_eq(1)
    w(2) = w(2)*baryc+(1._8-baryc) * w_eq(2)
    w(3) = w(3)*baryc+(1._8-baryc) * w_eq(3)

    w_check(1) = w(1)*baryc+(1._8-baryc) * w_eq(1)
    w_check(2) = w(2)*baryc+(1._8-baryc) * w_eq(2)
    w_check(3) = w(3)*baryc+(1._8-baryc) * w_eq(3)
       
    !if(maxval(abs(w-w_check))>1.e-14)then
    !  print *,'#warning in bgk_relax',maxval(abs(w-w_check))
    !  print *,'#w=',w
    !  print *,'#w_check=',w_check
    !endif   
      
  end subroutine bgk_relax
  
  subroutine compute_i0_and_alpha( &
    v, &
    dt, &
    xmin, &
    xmax, &
    N, &
    i0, &
    alpha)
    real*8, intent(in) :: v
    real*8, intent(in) :: dt
    real*8, intent(in) :: xmin
    real*8, intent(in) :: xmax
    integer, intent(in) :: N
    integer, intent(out) :: i0
    real*8, intent(out) :: alpha
    
    alpha = -v*dt/(xmax-xmin)
    alpha = alpha-floor(alpha)
    alpha = alpha*real(N,8)
    i0 = floor(alpha)
    if(i0 == N)then
      alpha = 0._8
      i0 = 0
    endif
    alpha = alpha-real(i0,8)
    
    
  end subroutine compute_i0_and_alpha

  subroutine compute_lag(lag,x,d)
    complex*16, intent(out) :: lag(:)
    complex*16, intent(in) :: x
    integer, intent(in) :: d
    
    integer :: i
    complex*16 :: a
    
    !buf(ix+i) <-> lag(j)
    !i=-d..d+1 <-> j=1..2*d+2 with j=i+d+1
    if(d>0)then 
    a=1._8
    do i=2,d
      a=a*(x*x-real(i,8)*real(i,8))/(real(d,8)*real(d,8))
    enddo
    a=a*(x+1._8)/real(d,8)
    a=a*(x-real(d,8)-1._8)/real(d,8)
    lag(d+1)=a*(x-1._8)/real(d,8)
    lag(d+2)=a*x/real(d,8)
    a=a*x*(x-1._8)/(real(d,8)*real(d,8))
    do i=-d,-1
      lag(i+d+1)=a/((x-real(i,8))/real(d,8))
    enddo  
    do i=2,d+1
      lag(i+d+1)=a/((x-real(i,8))/real(d,8))
    enddo
    a=1._8
    do i=-d,d+1 
      lag(i+d+1)=lag(i+d+1)*a
      a=a*real(d,8)/real(d+i+1,8)
    enddo
    a=1._8
    do i=d+1,-d,-1
      lag(i+d+1)=lag(i+d+1)*a
      a=a*real(d,8)/real(i-1-d-1,8)
    enddo
    else
      lag(1) = 1._8-x
      lag(2) = x
    endif
  
  end subroutine compute_lag

  subroutine semi_lag_advect(w,N,i0,lag,d)
    complex*16, intent(inout) :: w(:)
    integer, intent(in) :: N
    integer, intent(in) :: i0
    complex*16, intent(in) :: lag(:)
    integer, intent(in) :: d
    complex*16, allocatable :: buf(:)
    
    integer :: i
    integer :: j
    
    allocate(buf(N+1))
    
    buf(1:N+1) = w(1:N+1)
    
    !print *,real(w,8)
    
    do i=1,N+1
      w(i) = 0._8
      do j=-d,d+1
        w(i) = w(i)+lag(j+d+1)*buf(1+mod(i-1+j+i0+N,N))
      enddo
    enddo  
    !do i=1,N+1
    !  print *,i,real(buf(i),8),real(w(i),8)
    !enddo
    !stop  
  end subroutine semi_lag_advect
  
  !f_new(x) = f_old(x-I*dt)
  !method: f_new = fft(f_old)
  ! f_new(k) = f_new(k)*exp(-2*Pi*k*dt)
  !f(x) = exp(I*k*x*2*pi)
  !f(j/N) = exp(I*k*(j/N)*2*pi)
  !f(x-I*dt) = f(x)*exp(-k*2*pi*dt)
  subroutine semi_lag_advect_imag(w,N,dt)
    complex*16, intent(inout) :: w(:)
    integer, intent(in) :: N
    real*8, intent(in) :: dt
    
    type(C_PTR) :: fwd 
    type(C_PTR) :: bwd
    integer :: i
    integer :: Pi
    
    Pi = 3.1415926535897932385_8 
  
    !fwd = fftw_plan_dft_1d(N,w,w,FFTW_FORWARD,FFTW_ESTIMATE)
    !bwd = fftw_plan_dft_1d(N,w,w,FFTW_BACKWARD,FFTW_ESTIMATE)
    !call fftw_execute_dft(fwd, w, w)
    
    do i=1,N/2
      w(i) = w(i)*exp(-dt*2._8*Pi*real(i-1,8))
    enddo
    do i=N/2+1,N
      w(i) = w(i)*exp(-dt*2._8*Pi*real(i-1-N,8))
    enddo
    
    !call fftw_execute_dft(bwd, w, w)
    w = w/real(N,8)
    !call fftw_destroy_plan(fwd)
    !call fftw_destroy_plan(bwd)
    
  end subroutine semi_lag_advect_imag  

  subroutine compute_splitting_choice( &
    choice, &
    num_step, &
    begin_T, &
    dt_factor)
    character(len=*), intent(in) :: choice
    integer, intent(out) :: num_step
    logical, intent(out) :: begin_T
    complex*16, intent(out) :: dt_factor(:)
    select case (choice)
      case ("LIE_TV")
        begin_T = .true.
        num_step = 2
        if(size(dt_factor)<num_step)then
          print *,'#bad size for dt_factor'
          stop
        endif
        dt_factor(1) = 1._8
        dt_factor(2) = 1._8
      case ("COMPLEX2_TV")
        begin_T = .true.
        num_step = 4
        if(size(dt_factor)<num_step)then
          print *,'#bad size for dt_factor'
          stop
        endif
        dt_factor(1) = (0.5_8,0.5_8)
        dt_factor(2) = (0.5_8,0.5_8)
        dt_factor(3) = (0.5_8,-0.5_8)
        dt_factor(4) = (0.5_8,-0.5_8)
      case ("STRANG_TVT")
        begin_T = .true.
        num_step = 3
        if(size(dt_factor)<num_step)then
          print *,'#bad size for dt_factor'
          stop
        endif
        dt_factor(1) = 0.5_8
        dt_factor(2) = 1._8
        dt_factor(3) = 0.5_8
                
      case default
        print *,'#bad choice for splitting_choice'
        stop
    end select
  end subroutine  compute_splitting_choice

end module semi_lag_module
