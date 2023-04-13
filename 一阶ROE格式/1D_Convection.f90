module global
  implicit none                               ! 屏蔽掉Fortran中存在隐含的定义
   real::u(0:1000,0:1000)                     !  :: 表示修饰符的结束 定义二维数组u 数组大小1000*1000
end module
  
program main
  use global
  implicit none     
  
  ! 声明变量
  integer i,n,nt,nx                      ! 定义整数类型
  real k,dx,t,x,dt,Co
  character(len=10) method               ! 定义实数类型
  
    ! 输入数据
  k  = 1           
  x  = 2.0           ! y方向
  t  = 0.5           ! 时间(s)
  dx = 0.01          ! 网格大小
  Co = 0.5
  dt = Co*dx/k
  nt = int(t/dt)           ! 总的时间步数
  nx = int(x/dx)           ! 总的空间步数
  
! 初值及边界条件  
  do i = 0,nx
    if (i<=nx/2)then
      u(0,i) = 1.0
    else
      u(0,i) = 0.0
    end if      
  end do


  !选择插值方式
  write(*,*)"Please enter a difference method (LF  or    LW   or    BW   or    UW)?"
  read(*,*) method
  
  if(method == 'LF')then                                                                           ! LF
    do n = 0,nt-1                                                                                  ! 时间递进
      do i = 1,nx-1                                                                              ! x方向空间递进
        u(n+1,i) = 0.5*(u(n,i+1)+u(n,i-1)) - 0.5*Co*(u(n,i+1)-u(n,i-1))                            ! 代入差分方程  
      end do
      u(n+1,0)=u(n+1,1)
    end do

  else if(method == 'LW')then                                                                        ! LW
    do n = 0,nt-1                                                                                    ! 时间递进
      do i = 1,nx-1                                                                                  ! x方向空间递进
        u(n+1,i) = u(n,i) - Co*0.5*(u(n,i+1)-u(n,i-1)) + 0.5*Co**2*(u(n,i+1)-2*u(n,i)+u(n,i-1))      ! 代入差分方程
      end do
      u(n+1,0)=u(n+1,1)
    end do  

  else if(method == 'BW')then                                                                                  ! BW
    do n = 0,nt-1                                                                                              ! 时间递进
      do i = 2,nx                                                                                              ! x方向空间递进
        u(n+1,i) = u(n,i) - Co*0.5*(3*u(n,i)-4*u(n,i-1)+u(n,i-2)) + 0.5*Co**2*(u(n,i)-2*u(n,i-1)+u(n,i-2))     ! 代入差分方程
      end do
      u(n+1,1)=u(n+1,2)
      u(n+1,0)=u(n+1,1)
    end do  
    
  else if(method == 'UW')then                                                                                  ! BW
    do n = 0,nt-1                                                                                              ! 时间递进
      do i = 1,nx-1                                                                                              ! x方向空间递进
        u(n+1,i) = u(n,i) - Co*(u(n,i)-u(n,i-1))                                             ! 代入差分方程
      end do
      u(n+1,0)=u(n+1,1)
  end do   
  else
    write(*,*)"false"
    stop    
  end if

    ! 输出
  open(unit=11,file=trim(method)//'_result.txt')
  write(11,*)"    TIME           U            X"
  do n = 0,nt
    if(n*dt==0.0 .or. n==nt/4 .or.n==nt/2 .or.n==nt)then
      do i = 0,nx
        write(11,*)n*dt,u(n,i),i*dx-1     !读入这些时刻u的值。
      end do
      write(11,*)
    end if
  end do

  close(11)
  

  stop 
  end program    