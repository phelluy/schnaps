!gfortran -O3 -I/usr/local/include /usr/local/lib/libfftw3.dylib semi_lag_module.F90 semi_lag.F90 -o semi_lag
!gfortran -O3 semi_lag_module.F90 test_diagnostics_module.F90 semi_lag.F90 -o semi_lag
program semi_lag
  use semi_lag_module
  use test_diagnostics_module
 
implicit none

  real*8 :: xmin
  real*8 :: xmax
  real*8 :: cfl
  complex*16 :: c
  real*8 :: tmax
  integer :: M
  integer :: N
  integer :: d
  
  !sound speed
  real*8 :: CSON
  !estimated maximal wave speed
  real*8 :: VMAX
  
  real*8 :: t
  complex*16 :: dt
  complex*16 :: z
  real*8 :: dx
  complex*16, allocatable :: w(:,:)
  integer :: i
  real *8 :: x
  integer :: file_id
  real*8, allocatable :: v(:)
  real*8, allocatable :: alpha(:)
  integer, allocatable :: i0(:)
  complex*16, allocatable :: lag(:,:)
  integer :: num
  logical :: smooth
  complex*16 :: baryc
  real*8 :: tau
  real*8 :: relax
  integer :: num_step
  real*8 :: dt_big
  integer :: sub_step
  integer :: num_sub_step
  complex*16, allocatable :: dt_factor(:)
  integer :: num_sub_step_max
  logical :: begin_T
  logical :: split_T
  character(len=256) :: splitting_choice
  real*8 :: wave_factor
  real*8 :: factor_pos,factor_neg

  
   
  smooth = .true.
  file_id = 99
  num_sub_step_max = 100
    
  xmin = -2._8
  xmax = 2._8
  cfl = 0.1_8
  c = (0.5_8,0.5_8)
  !c = (0.5_8,0._8)
  M = 3
  N= 1000
  d = 8
  tau =0.0e-12_8

  wave_factor = 1.0_8
  tmax = 0.1_8
  
  VMAX = 2.0_8 
  relax = 1._8/tau
  splitting_choice = "COMPLEX2_TV"
  splitting_choice = "STRANG_TVT"
  
  print *,'#xmin=',xmin
  print *,'#xmax=',xmax
  print *,'#cfl=',cfl
  print *,'#c=',c
  print *,'#tmax=',tmax
  print *,'#M=',M
  print *,'#CSON=',CSON
  print *,'#VMAX=',VMAX
  
  t = 0.0_8
  dx = (xmax-xmin)/real(N,8)
  
  allocate(w(M,N+1))
  allocate(v(M))
  allocate(alpha(M))
  allocate(i0(M))
  allocate(lag(2*d+2,M))
  allocate(dt_factor(num_sub_step_max))
  v(1) = -VMAX
  v(2) = 0._8
  v(3) = VMAX 
  
  dt = 0.005_8
  num_step = floor(tmax/real(dt,8)+1.e-10_8)
  
  dt_big = dt
  !call compute_dt( &
  !  VMAX, &
  !  cfl, &
  !  dx, &
  !  c, &
  !  dt)
  !call compute_dt_big( &
  !  VMAX, &
  !  cfl, &
  !  dx, &
  !  dt_big)

  t = 0._8
  if(wave_factor .eq. 0.0_8)then
     do i=1,N+1
        x = xmin+real(i-1,8)*dx
        if(smooth)then
           CSON = 1._8
           call init_sol_Wave_smooth( &
                x, &
                t, &
                w(:,i), &
                VMAX)
        else
           CSON = 1._8
           call init_sol_Wave_nonsmooth( &
                x, &
                t, &
                w(:,i), &
                VMAX)         
        endif
     enddo
  else
     do i=1,N+1
        x = xmin+real(i-1,8)*dx
        if(smooth)then
            CSON = 0.6_8
            call init_sol_Euler_smooth( &
                x, &
                t, &
                w(:,i), &
                VMAX=VMAX, &
                CSON=CSON)
         else
            CSON = 0.6_8
            call init_sol_Euler_nonsmooth( &
                 x, &
                 t, &
                 w(:,i), &
                 VMAX=VMAX, &
                 CSON=CSON )
         endif
      end do
   endif
     
  baryc = 0._8
  do i=1,N+1
    call bgk_relax(w(:,i),baryc,VMAX, CSON, wave_factor)
  enddo
  
  num=0

  
  
  call compute_splitting_choice( &
    splitting_choice, &
    num_sub_step, &
    begin_T, &
    dt_factor)
    
  num_sub_step = 2
  begin_T = .true.
  dt_factor(1) = 1._8
  dt_factor(2) = 1._8
  
  !num_sub_step = 3
  !begin_T = .true.
  !dt_factor(1) = 0.5_8
  !dt_factor(2) = 1._8
  !dt_factor(3) = 0.5_8

  !num_sub_step = 7
  !begin_T = .true.
  !dt_factor(1) = 0.675603595979829_8
  !dt_factor(2) = 1.351207191959658_8
  !dt_factor(3) = -0.17560359597982855_8
  !dt_factor(4) = -1.702414383919315_8
  !dt_factor(5) = dt_factor(3)
  !dt_factor(6) = dt_factor(2)
  !dt_factor(7) = dt_factor(1)

  !num_sub_step = 11
  !begin_T = .true.
  !factor_pos = 0.41449077179_8	
  !factor_neg = 0.65796308717_8
  !dt_factor(1) = 0.5_8 * factor_pos
  !dt_factor(2) = factor_pos
  !dt_factor(3) = factor_pos
  !dt_factor(4) = factor_pos
  !dt_factor(5) = 0.5_8 *factor_pos - 0.5_8 * factor_neg
  !dt_factor(6) = factor_neg
  !dt_factor(7) = 0.5_8 *factor_pos - 0.5_8 * factor_neg
  !dt_factor(8) = factor_pos
  !dt_factor(9) = factor_pos
  !dt_factor(10) = factor_pos
  !dt_factor(11) = 0.5 * factor_pos


  !dt_factor(3) = conjg(c)
  !dt_factor(4) = conjg(c)
  
  print *,'#num_step=',num_step
  
  !do while (t<tmax-1e-10_8)
  do num=1,num_step
    split_T = begin_T
    do sub_step=1,num_sub_step
      dt = dt_factor(sub_step)*dt_big
      if(split_T)then
        !semi_lag step
        do i=1,M
          call compute_i0_and_alpha( &
            v(i), &
            real(dt,8), &
            xmin, &
            xmax, &
            N, &
            i0(i), &
            alpha(i))
        enddo
        if(t==0._8)then
          print *,'#i0=',i0
          print *,'#alpha=',alpha
        endif
        do i=1,M
          call compute_lag(lag(:,i),alpha(i)+(0._8,0._8),d)
          !print *,'#lag=',lag(:,i)
        enddo  
        do i=1,M
          call semi_lag_advect(w(i,:),N,i0(i),lag(:,i),d)
        enddo
        !do i=1,M
        !  call semi_lag_advect_imag(w(i,:),N,v(i)*imag(dt))
        !enddo  
      else
        do i=1,M
          !call semi_lag_advect_imag(w(i,:),N,v(i)*imag(dt))
        enddo  
        !relaxation step
        if(tau == 0) then
           baryc =-1._8
        else
           z= tau + 0.5_8 * dt
           baryc = (tau - 0.5_8 * dt) / z
        end if
        if(num_sub_step == 2) then
           if(tau == 0) then
              baryc =0.0_8
           else
              z= tau + dt
              baryc = (tau - dt) / z
           end if
        end if
        !baryc = 0._8
        !baryc = 0.5_8*dt
        !baryc = 1._8/(1._8+relax*dt)    
        !print *,'#baryc=',baryc
        do i=1,N+1
           call bgk_relax(w(:,i),baryc,VMAX, CSON, wave_factor)
        enddo
        !do i=1,M
          !call semi_lag_advect_imag(w(i,:),N,-v(i)*imag(dt))
        !enddo  
      endif
      split_T = .not.(split_T)
    enddo        
    t = t+dt_big
    print *,'#t=',t,num
  enddo

  baryc = 0._8
  do i=1,N+1
    call bgk_relax(w(:,i),baryc,VMAX, CSON, wave_factor)
  enddo


  if(wave_factor == 0.0_8)then
     open(unit = file_id,file='solex.dat')
     call write_wave_solex( &
          file_id, &
          xmin, &
          xmax, &
          w, &
          VMAX, &
          smooth, &
          tmax)
  endif
  
  open(unit = file_id,file='semi_lag.dat')
  call write_sol( &
          file_id, &
          xmin, &
          xmax, &
          VMAX, &
          w)

  
end program semi_lag
