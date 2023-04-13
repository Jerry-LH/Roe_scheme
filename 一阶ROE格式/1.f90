program test
    implicit none
    integer::i
    real :: a, b, c, foo ,d,min_mod,dx
    do i=0,10
    b = -5.0
    c=-3.0
    a = min_mod(b+3.0,c)
    print *,a,b,c
    enddo
    
    dx = min_mod(b, c)
    print*,dx
end program

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