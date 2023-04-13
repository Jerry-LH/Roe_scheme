module global
  implicit none
  integer::i,j,m,l,rk,x_step=200
  real*8 ::x(0:200),t=0,dt=1e-5,gama=1.4,U(0:201,3),lamda(3,3)=0,dx,&
        u_a,H_a,c_a,k(3,3),k_inverse(3,3),F_half(0:201,3),vis(3,3),U00(0:201,3),U01(0:201,3),U02(0:201,3),vis1                      
end module

program main
  use global
  do i=0,x_step
    dx=1./x_step                           ! 网格划分
    x(i)=-0.5*dx+dx*i
    if(x(i)<=0.3)then                      ! 初始条件
      U(i,1)=1
      U(i,2)=0
      U(i,3)=2.5
    else
      U(i,1)=0.125
      U(i,2)=0
      U(i,3)=0.25
    end if
  end do

  do i=1,60000
    U00=U
    do rk=1,3
      call bound                                  ! 边界无反射条件
      call flux
      call update
    end do
    t=t+dt
  end do
  call output
  stop
end program

subroutine bound
  use global
  do j=0,x_step
    F_half(j,1)=U(j,2)
    F_half(j,2)=U(j,2)**2/U(j,1)+(gama-1)*(U(j,3)-0.5*U(j,2)**2/U(j,1))
    F_half(j,3)=(U(j,3)+(gama-1)*(U(j,3)-0.5*U(j,2)**2/U(j,1)))*U(j,2)/U(j,1)
  end do
  U(0,:)=U(1,:)
  U(x_step+1,:)=U(x_step,:)
  F_half(0,:)=F_half(1,:)
  F_half(x_step+1,:)=F_half(x_step,:)
end subroutine bound

subroutine flux
  use global
  do j=0,x_step

    u_a=(sqrt(U(j,1))*U(j,2)/U(j,1)+sqrt(U(j+1,1))*U(j+1,2)/U(j+1,1))/(sqrt(U(j,1))+sqrt(U(j+1,1)))
    H_a=(sqrt(U(j,1))*((U(j,3)+F_half(j,2)-0.5*U(j,1)*(U(j,2)/U(j,1))**2))/U(j,1)+sqrt(U(j+1,1))*&
        ((U(j+1,3)+F_half(j+1,2)-0.5*U(j+1,1)*(U(j+1,2)/U(j+1,1))**2))/U(j+1,1))/(sqrt(U(j,1))+sqrt(U(j+1,1)))
    c_a=sqrt((gama-1)*(H_a-0.5*u_a**2))

    lamda(1,1)=abs(u_a-c_a)
    lamda(2,2)=abs(u_a)
    lamda(3,3)=abs(u_a+c_a)

    k(1,1)=1                                                                                ! 特征向量矩阵定义
    k(1,2)=1
    k(1,3)=1
    k(2,1)=u_a-c_a
    k(2,2)=u_a
    k(2,3)=u_a+c_a
    k(3,1)=H_a-u_a*c_a
    k(3,2)=0.5*u_a**2
    k(3,3)=H_a+u_a*c_a

    k_inverse(1,1)=(gama-1)*u_a**2/(4*c_a**2)+0.5*u_a/c_a                   ! 特征向量逆矩阵定义
    k_inverse(1,2)=-(gama-1)*u_a/(2*c_a**2)-1./(2*c_a)
    k_inverse(1,3)=(gama-1)/(2*c_a**2)
    k_inverse(2,1)=1-(gama-1)*u_a**2/(2*c_a**2)
    k_inverse(2,2)=(gama-1)*u_a/(c_a**2)
    k_inverse(2,3)=-(gama-1)/(c_a**2)
    k_inverse(3,1)=(gama-1)*u_a**2/(4*c_a**2)-0.5*u_a/c_a
    k_inverse(3,2)=-(gama-1)*u_a/(2*c_a**2)+1./(2*c_a)
    k_inverse(3,3)=(gama-1)/(2*c_a**2)
    
    vis=MATMUL(MATMUL(k_inverse,lamda),k)
    
    do m=1,3
      vis1=0
      do l=1,3
        vis1=vis1+vis(l,m)*(U(j+1,l)-U(j,l))
      end do
      F_half(j,m)=0.5*(F_half(j,m)+F_half(j+1,m))-0.5*vis1
    end do

  end do
end subroutine flux

subroutine update    
  use global
  if(rk.eq.1)then
    do j=1,x_step                                                                            ! 时间推进RK3法
      U01(j,:)=U00(j,:)-dt/dx*(F_half(j,:)-F_half(j-1,:))
    end do
    U=U01

  else if(rk.eq.2)then
    do j=1,x_step      
      U02(j,:)=3./4*U00(j,:)+1./4*U01(j,:)-1./4*dt/dx*(F_half(j,:)-F_half(j-1,:))
    end do
    U=U02

  else if(rk.eq.3)then
    do j=1,x_step
      U(j,:)=1./3*U00(j,:)+2./3*U02(j,:)-2./3*dt/dx*(F_half(j,:)-F_half(j-1,:))
    end do

  end if

end subroutine update

subroutine output
  use global
  open(1,file='out.dat')
  write(1,*)'Variables="X","Density","Velocity","Pressure","Energy"'

  do i=1,x_step
    write(1,*)x(i),U(i,1),U(i,2)/U(i,1),(gama-1)*(U(j,3)-0.5*U(i,2)**2/U(i,1)),U(i,3)/U(i,1)
  end do

end subroutine output