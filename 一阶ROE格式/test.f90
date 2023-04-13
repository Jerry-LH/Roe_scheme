module global
  implicit none
  integer::i,j,m,l,rk,x_step=200
  real::x(0:200),t=0,dt=1e-5,gama=1.4,U(-1:202,3),F(0:201,3),lamda(3),dx,&
        u_a,h_a,c_a,k(3,3),F_half(0:201,3),U00(-1:202,3),U01(-1:202,3),U02(-1:202,3),alpha(3)                    
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

  do i=1,20000
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
  U(0,:)=U(1,:)
  U(-1,:)=U(0,:)
  U(x_step+1,:)=U(x_step,:)
  U(x_step+2,:)=U(x_step+1,:)
end subroutine bound

subroutine flux
  use global
  real::rl,rr,ul,ur,el,er,pl,pr,hl,hr,pll,prr,pl0,pr0,min_mod
  do j=0,x_step

    rl=U(j,1)
    rr=U(j+1,1)
    ul=U(j,2)/U(j,1)
    ur=U(j+1,2)/U(j+1,1)
    el=U(j,3)/U(j,1)
    er=U(j+1,3)/U(j+1,1)
    pl0=(gama-1)*(rl*el-0.5*rl*ul**2)
    pr0=(gama-1)*(rr*er-0.5*rr*ur**2)
    pll=(gama-1)*(U(j-1,1)*U(j-1,3)/U(j-1,1)-0.5*U(j-1,1)*U(j-1,2)/U(j-1,1)**2)
    prr=(gama-1)*(U(j+2,1)*U(j+2,3)/U(j+2,1)-0.5*U(j+2,1)*U(j+2,2)/U(j+2,1)**2)
    
    rl=rl+0.5*min_mod(U(j+1,1)-U(j,1),U(j,1)-U(j-1,1))
    rr=rr-0.5*min_mod(U(j+1,1)-U(j,1),U(j+2,1)-U(j+1,1))
    ul=ul+0.5*min_mod(U(j+1,2)/U(j+1,1)-U(j,2)/U(j,1),U(j,2)/U(j,1)-U(j-1,2)/U(j-1,1))
    ur=ur-0.5*min_mod(U(j+1,2)/U(j+1,1)-U(j,2)/U(j,1),U(j+2,2)/U(j+2,1)-U(j+1,2)/U(j+1,1))
    pl=pl0+0.5*min_mod(pr0-pl0,pl0-pll)
    pr=pr0-0.5*min_mod(pr0-pl0,prr-pr0)
    el=(pl/(gama-1)+0.5*rl*ul**2)/rl
    er=(pr/(gama-1)+0.5*rr*ur**2)/rr

    hl=el+pl/rl
    hr=er+pr/rr

    F_half(j,1)=0.5*(rl*ul+rr*ur)
    F_half(j,2)=0.5*(rl*ul**2+pl+rr*ur**2+pr)
    F_half(j,3)=0.5*(rl*hl*ul+rr*hr*ur)

    u_a=(sqrt(rl)*ul+sqrt(rr)*ur)/(sqrt(rr)+sqrt(rl))                                        ! 取该点的左右两值，方便下面计算u_j+1/2
    h_a=(sqrt(rl)*hl+sqrt(rr)*hr)/(sqrt(rr)+sqrt(rl))
    c_a=sqrt((gama-1)*(h_a-0.5*u_a**2))


    lamda(1)=abs(u_a-c_a)
    lamda(2)=abs(u_a)
    lamda(3)=abs(u_a+c_a)

    k(1,1)=1                                                                                ! 特征向量矩阵定义
    k(1,2)=1
    k(1,3)=1
    k(2,1)=u_a-c_a
    k(2,2)=u_a
    k(2,3)=u_a+c_a
    k(3,1)=h_a-u_a*c_a
    k(3,2)=0.5*u_a**2
    k(3,3)=h_a+u_a*c_a

    alpha(2)=(gama-1)/(c_a**2)*((rr-rl)*(h_a-u_a**2)+u_a*(rr*ur-rl*ul)-(rr*er-rl*el))
    alpha(1)=0.5*((rr-rl)*(u_a+c_a)-(rr*ur-rl*ul)-c_a*alpha(2))/c_a
    alpha(3)=(rr-rl)-(alpha(1)+alpha(2))

    F_half(j,1)=F_half(j,1)-0.5*(lamda(1)*alpha(1)*k(1,1)+lamda(2)*alpha(2)*k(1,2)+lamda(3)*alpha(3)*k(1,3))
    F_half(j,2)=F_half(j,2)-0.5*(lamda(1)*alpha(1)*k(2,1)+lamda(2)*alpha(2)*k(2,2)+lamda(3)*alpha(3)*k(2,3))
    F_half(j,3)=F_half(j,3)-0.5*(lamda(1)*alpha(1)*k(3,1)+lamda(2)*alpha(2)*k(3,2)+lamda(3)*alpha(3)*k(3,3))
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
    write(1,*)x(i),U(i,1),U(i,2)/U(i,1),(gama-1)*(U(i,3)-0.5*U(i,1)*(U(i,2)/U(i,1))**2),U(i,3)/U(i,1)
  end do

end subroutine output

function min_mod(input1,input2)
    implicit none
    real :: min_mod
    real, intent(in) :: input1
    real, intent(out) :: input2
    if(input1*input2<=0)then
      min_mod=0
    else if(input1*input2>0)then
      if(abs(input1)>=abs(input2))then
        min_mod=input2
      else
        min_mod=input1
      end if
    end if
end function