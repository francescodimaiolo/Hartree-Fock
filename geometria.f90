program geometria
  implicit none
  real*8:: l, l1, alpha
  integer:: i
  real*8:: x, y, pi
  l1=1.50
  l=l1
  pi=dacos(-1.d0)
  alpha=pi/6
  write(*,*) dsin(alpha), dcos(alpha)
  open(1,file='geometrie.dat')
  do i=1,13
     if(i==1)then
        x=-l*dcos(alpha)
        y=l1+l*dsin(alpha)
     endif
     if(i==2)then
        x=-2*l*dcos(alpha)
        y=l
     endif
     if(i==3)then
        x=-(l+l1)*dcos(alpha)
        y=0
     endif
     if(i==4)then
        x=-l1*dcos(alpha)
        y=-l1*dsin(alpha)
     endif
     if(i==5)then
        x=-l*dcos(alpha)
        y=-(l1*dsin(alpha)+l)
     endif
     if(i==6)then
        x=0
        y=-(l1*dsin(alpha)+l+l*dsin(alpha))
     endif
     if(i==7)then
        x=l*dcos(alpha)
        y=-(l1*dsin(alpha)+l)
     endif
     if(i==8)then
        x=l1*dcos(alpha)
        y=-(l1*dsin(alpha))
     endif
     if(i==9)then
        x=(l+l1)*dcos(alpha)
        y=0
     endif
     if(i==10)then
        x=2*l*dcos(alpha)
        y=l
     endif
     if(i==11)then
        x=(l)*dcos(alpha)
        y=(l1+l*dsin(alpha))
     endif
     if(i==12)then
        x=0
        y=(l1)
     endif
     if(i==13)then
        x=0
        y=0
     endif
     write(1,*) x, y, i
  enddo
  close(1)
end program geometria
