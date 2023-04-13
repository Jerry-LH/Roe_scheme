module global
  implicit none       ! 屏蔽掉Fortran中存在隐含的定义
  real::u(1000,1000)  !  :: 表示修饰符的结束 定义二维数组u 数组大小1000*1000
end module
  
program main
  use global
  implicit none     
  
  ! 声明变量
  integer i,n,h,nt,ny               ! 定义整数类型
  real nu,dt,dy,t,m                 ! 定义实数类型
  character(len=20) chazhi          ! 定义chazhi为内在数据类型
  real::e(1000),b(1000),c(1000)     ! 定义一维数组e,b,c 数组大小均为1000
  

  ! 输入数据
  nu = 0.000217      ! 扩散系数 单位m**2/s
  ny = 41            ! y方向
  t  = 1.08          ! 时间(s)
  dy = 0.001         ! 网格大小
  dt = 0.002         ! 时间步长
  nt = int(t/dt)     ! 总的时间步数
  m  = nu*dt/dy**2   ! 定义m为差分方法中的系数
   
  ! 初值及边界条件
  do i = 1,ny        ! 做循环，使矩阵u的第1列，第1至41项初始为0
    u(i,1) = 0.0
  end do

  do n = 1,nt        ! 定义y=0时，u=U0 y=h时，u=0 这里t=0，y=0的初始条件重合了，所以上面不用在重新设置
    u(1,n)  = 40.0
    u(ny,n) = 0.0
  end do
  
  !选择插值方式
  write(*,*)"FTCS  or    DF   or    LA   or   CN  ?"
  read(*,*) chazhi
  
  if(chazhi == 'FTCS')then                                       ! FTCS
    do n = 1,nt-1                                                ! 时间递进
      do i = 2,ny-1                                              ! y方向空间递进
        u(i,n+1) = u(i,n) + m *(u(i+1,n)-2*u(i,n)+u(i-1,n))      ! 代入差分方程
      end do
    end do
    
  else if(chazhi == 'DF')then  ! DF
    do n = 2,nt-1  
      do i = 2,ny-1 
        u(i,n+1) = (1-2*m)/(1+2*m) * u(i,n-1) + 2*m/(1+2*m) *(u(i+1,n)+u(i-1,n))
      end do
    end do
    
  else if(chazhi == 'LA')then  ! LA
      do i = 1, 1000
        e(i)  = m              ! LA方法系数A,i-1
        b(i)  = -(2*m+1)       ! LA方法系数B,i
        c(i)  = m              ! LA方法系数C,i+1
      end do
      
      do n = 1,nt-1
        call LAsolve(e,b,c,ny,n)  ! 求解n+1时刻的三对角方程组
      end do
      
  else if(chazhi == 'CN')then  ! CN
      do i = 1, 1000
        e(:)  = m/2            ! LA方法系数A,i-1
        b(:)  = -(m+1)         ! LA方法系数B,i
        c(:)  = m/2            ! LA方法系数C,i+1
      end do
      
      do n = 1,nt-1
        call CNsolve(e,b,c,ny,n,m) ! 求解n+1时刻的三对角方程组
      end do
  else
    write(*,*)"false"
    stop
    
  end if
  
  ! 输出
  open(unit=100,file=trim(chazhi)//'_result.txt')
  write(100,*)"    TIME           U            Y"
  do n = 1,nt
    if(n*dt==0.18 .or. n*dt==0.36.or.n*dt==0.54.or.n*dt==0.72.or.n*dt==0.90.or.n*dt==1.08)then
      do i = 1,ny
        write(100,*)n*dt,u(i,n),(i-1)*dy!读入这些时刻u的值。
      end do
      write(100,*)
    end if
  end do

  close(100)
  
  stop 
  end program

subroutine  LAsolve(e,b,c,ny,n)
  use global
  implicit none
  
  integer i,k
  integer n,ny
  real::e(1000),b(1000),c(1000),y(1000),f(1000)
  real::L(1000),M(1000)
    
  M(1) = b(1)
  do i=2,ny
    L(i)=e(i)/M(i-1)
    M(i)=b(i)-L(i)*c(i-1)
  end do
  f(1) = -u(1,n)-e(1)*u(1,n)        ! 1处的边界值
  f(ny) = -u(ny,n)-c(ny)*u(ny,n)    ! ny处的边界值
  do i = 2,ny-1
    f(i) = -u(i,n)
  end do

    
  !------开始回带,求得y,1>>ny
  y(1)=f(1)
  do i=2,ny
    y(i)=f(i)-L(i)*y(i-1)
  end do
  
  !-----开始回带，求得n+1时刻的u,ny>>1
  u(ny-1,n+1)=y(ny)/m(ny)   
  do i=ny-1,2,-1
    u(i,n+1)=(y(i)-c(i)*u(i+1,n))/m(i)
  end do
end subroutine LAsolve
  
  
subroutine  CNsolve(e,b,c,ny,n,w)
  use global
  implicit none
  
  integer i,k
  integer n,ny
  real w
  real::e(1000),b(1000),c(1000),y(1000),f(1000)
  real::L(1000),M(1000)
    
  M(2) = b(2)
  do i=3,ny-1
    L(i)=e(i)/M(i-1)
    M(i)=b(i)-L(i)*c(i-1)
  end do
                                                                       
  f(2)    = (w-1)*u(2,n) - w/2*(u(3,n)+u(1,n))-e(1)*u(1,n)           ! 1处的边界值
  f(ny-1) = (w-1)*u(ny-1,n) - w/2*(u(ny,n)+u(ny-1,n))-c(ny)*u(ny,n)  ! ny处的边界值
  do i = 3,ny-2
    f(i) = (w-1)*u(i,n) - w/2*(u(i+1,n)+u(i-1,n))
  end do
    
  !------开始回带,求得y,1>>ny
  y(2)=f(2)
  do i=3,ny-1
    y(i)=f(i)-L(i)*y(i-1)
  end do
  
  !-----开始回带，求得n+1时刻的u,ny>>1
  u(ny-1,n+1)=y(ny-1)/m(ny-1)   
  do i=ny-2,2,-1
    u(i,n+1)=(y(i)-c(i)*u(i+1,n))/m(i)
  end do
  
end subroutine CNsolve