module global
  implicit none                ! 屏蔽掉Fortran中存在隐含的定义
   real*8::u(0:10000,0:10000)     !  :: 表示修饰符的结束 定义二维数组u 数组大小1000*1000
end module
  
program main
  use global
  implicit none     
  
  ! 声明变量
  integer i,n,xc,tc                ! 定义整数类型
  real*8 k,dx,t,x,m,PI,nt,nx,dt                ! 定义实数类型
  
    ! 输入数据
  PI = DACOS(-1.D0)  ! PI
  k = 0.5           ! 扩散系数 单位m**2/s
  x  = 1.0           ! y方向
  t  = 0.2           ! 时间(s)
  dx = 0.01          ! 网格大小

  write(*,*)"Please enter a time step"
  read(*,*) dt
 
  nt = t/dt          ! 总的时间步数
  tc = int(nt)
  nx = x/dx           ! 总的空间步数
  xc = int(nx)
  m  = k*dt/dx**2   ! 定义m为差分方法中的系数
  
    ! 初值及边界条件
  do i = 0,xc
    u(0,i) = i/nx - 1 + sin(4*PI*(i/nx))
  end do

  do n = 0,tc
    u(n,0) = -1
    u(n,xc) = 0
  end do
  
  do n = 0,tc-1                                                ! 时间递进
    do i = 1,xc-1                                            ! x方向空间递进
      u(n+1,i) = u(n,i) + m*(u(n,i+1)-2*u(n,i)+u(n,i-1))     ! 代入差分方程
    end do
  end do
  
    ! 输出
  open(unit=11,file='1_result.txt')
  write(11,*)"    TIME           U            X"
  do n = 0,tc
    do i = 0,xc
      write(11,*)n*dt,u(n,i),i*dx     !读入这些时刻u的值。
    end do
    write(11,*)
  end do

  close(11)
  

  stop 
  end program    