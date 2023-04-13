module gx   
    integer,parameter::m=200                                                                ! 定义空间步数、比热容比、时间步长
    real,parameter::gama=1.4,dt=1e-5
end module gx
    
module gx1                                                                                  ! 算出空间步长，设置初始时间
    use gx
    real::dx=1./m    
    real::t=0      
    integer rk                                                                              ! 设置循环系数
    real::x(0:m+1)                                                                          ! 设置空间步长矩阵
    real::r(0:m+1),ru(0:m+1),re(0:m+1),r00(0:m+1),r01(0:m+1),r02(0:m+1),&                   ! r 密度    ru 密度*速度    re 
        ru00(0:m+1),ru01(0:m+1),ru02(0:m+1),re00(0:m+1),re01(0:m+1),re02(0:m+1)
    real::rhs(0:m,1:3)
    real::p(0:m),u(0:m),e(0:m)
end module gx1
      
program yiwei    
    use gx
    use gx1
    integer::i,j
    call mesh
    call initial

    do j=1,15000                                                                              ! 时间推进，0.15s
        do i=1,m
            r00(i)=r(i)
            ru00(i)=ru(i)
            re00(i)=re(i)
        end do
        


        do rk=1,3
            call bound                                                                                 ! 边界无反射条件 
            call flux 
            call update
        end do
        
        t=t+dt

    end do

    call output
    
    stop
    
end program yiwei

subroutine mesh    
    use gx
    use gx1

    do i=0,m+1
        x(i)=dx*0.5+(i-1)*dx
    end do

    return
end subroutine mesh
      
subroutine initial
    use gx
    use gx1
    do i=1,m
        if(x(i)<=0.3) then
            r(i)=1.0
            u(i)=0.
            p(i)=1.0
            ru(i)=r(i)*u(i)
            re(i)=p(i)/(gama-1)+0.5*r(i)*ru(i)**2                                  ! 补充方程
        else
            r(i)=0.125
            u(i)=0.
            p(i)=0.1
            ru(i)=r(i)*u(i)
            re(i)=p(i)/(gama-1)+0.5*r(i)*ru(i)**2
        end if
    end do
    return
end subroutine initial

subroutine bound     
    use gx
    use gx1
    real::rl,rr,rul,rur,rel,rer,ul,ur,pl,pr,hl,hr,dr,dru,dre,ua,ha,ca
        
    r(0)=r(1)
    ru(0)=ru(1)
    re(0)=re(1)
    r(m+1)=r(m)
    ru(m+1)=ru(m)
    re(m+1)=re(m)
    
    return
end subroutine bound
      
subroutine flux 
    use gx
    use gx1
    real::k(3,3),lamda(3),alpha(3),w(3),W1(3),beta(3),r0(3)
    real::rl,rr,rul,rur,rel,rer,ul,ur,pl,pr,hl,hr,dr,dru,dre,ua,ha,ca
    do i=0,m
        rl=r(i)                                                                               ! 定义左右侧的守恒变量
        rr=r(i+1)
        rul=ru(i)
        rur=ru(i+1)
        rel=re(i)
        rer=re(i+1)

        ul=rul/rl
        ur=rur/rr
        pl=(gama-1)*(rel-0.5*rl*ul**2)
        pr=(gama-1)*(rer-0.5*rr*ur**2)
        hl=(rel+pl)/rl
        hr=(rer+pr)/rr
     

        rhs(i,1)=0.5*(rul+rur)                                                                 ! 中心性有限体积法
        rhs(i,2)=0.5*(rul*ul+pl+rur*ur+pr)
        rhs(i,3)=0.5*(rl*hl*ul+rr*hr*ur)
        
        dr=rr-rl                                                                               ! 差值delta（F）
        dru=rur-rul
        dre=rer-rel                                                                            
        
        ua=(sqrt(rl)*ul+sqrt(rr)*ur)/(sqrt(rr)+sqrt(rl))                                        ! 取该点的左右两值，方便下面计算u_j+1/2
        ha=(sqrt(rl)*hl+sqrt(rr)*hr)/(sqrt(rr)+sqrt(rl))
        ca=sqrt((gama-1)*(ha-0.5*ua**2))
        
        lamda(1)=abs(ua-ca)                                                                     ! roe特征值取绝对值 
        lamda(2)=abs(ua)
        lamda(3)=abs(ua+ca)
    
        k(1,1)=1                                                                                ! 特征向量矩阵定义
        k(1,2)=ua-ca
        k(1,3)=ha-ua*ca
        k(2,1)=1
        k(2,2)=ua+ca
        k(2,3)=0.5*ua**2
        k(3,1)=1
        k(3,2)=ua+ca
        k(3,3)=ha+ua*ca
    
        alpha(2)=(gama-1)/(ca**2)*(dr*(ha-ua**2)+ua*dru-dre)
        alpha(1)=0.5*(dr*(ua+ca)-dru-ca*alpha(2))/ca
        alpha(3)=dr-(alpha(1)+alpha(2))
    
        rhs(i,1)=rhs(i,1)-0.5*(lamda(1)*alpha(1)*k(1,1)+lamda(2)*alpha(2)*k(2,1)+lamda(3)*alpha(3)*k(3,1))          ! 通过roe格式求出半个空间步长的值
        rhs(i,2)=rhs(i,2)-0.5*(lamda(1)*alpha(1)*k(1,2)+lamda(2)*alpha(2)*k(2,2)+lamda(3)*alpha(3)*k(3,2))
        rhs(i,3)=rhs(i,3)-0.5*(lamda(1)*alpha(1)*k(1,3)+lamda(2)*alpha(2)*k(2,3)+lamda(3)*alpha(3)*k(3,3))
    end do
    return
end subroutine flux

subroutine update    
    use gx
    use gx1

    if(rk.eq.1) then
        do i=1,m                                                                            ! 时间推进RK3法
            r01(i)=r00(i)-dt/dx*(rhs(i,1)-rhs(i-1,1))
            ru01(i)=ru00(i)-dt/dx*(rhs(i,2)-rhs(i-1,2))
            re01(i)=re00(i)-dt/dx*(rhs(i,3)-rhs(i-1,3))
            r(i)=r01(i)
            ru(i)=ru01(i)
            re(i)=re01(i)
        end do

    else if(rk.eq.2)  then
        do i=1,m
            r02(i)=3./4*r00(i)+1./4*r01(i)-1./4*dt/dx*(rhs(i,1)-rhs(i-1,1))
            ru02(i)=3./4*ru00(i)+1./4*ru01(i)-1./4*dt/dx*(rhs(i,2)-rhs(i-1,2))
            re02(i)=3./4*re00(i)+1./4*re01(i)-1./4*dt/dx*(rhs(i,3)-rhs(i-1,3))
            r(i)=r02(i)
            ru(i)=ru02(i)
            re(i)=re02(i)
        end do

    else if(rk.eq.3) then
        do i=1,m
            r(i)=1./3*r00(i)+2./3*r02(i)-2./3*dt/dx*(rhs(i,1)-rhs(i-1,1))
            ru(i)=1./3*ru00(i)+2./3*ru02(i)-2./3*dt/dx*(rhs(i,2)-rhs(i-1,2))
            re(i)=1./3*re00(i)+2./3*re02(i)-2./3*dt/dx*(rhs(i,3)-rhs(i-1,3))
        end do

    end if

end subroutine update
    
subroutine output
    use gx
    use gx1
    open(1,file='out1.dat')
    write(1,*)'Variables="X","Density","Velocity","Pressure","Energy"'

    do i=1,m
        p(i)=(gama-1)*(re(i)-0.5*ru(i)**2/r(i))
        u(i)=ru(i)/r(i)
        e(i)=re(i)/r(i)-0.5*u(i)**2
        write(1,*)x(i),r(i),u(i),p(i),e(i)
    end do

end subroutine output
      



























      