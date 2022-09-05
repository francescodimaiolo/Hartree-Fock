!!!! trovare configurazione massima e minima in bit representation
program basis
  implicit none
  integer:: nso,i,max,min, n, max1,max2, min1, min2, e
  double precision::o
  write(*,*) 'dammi il numero elettroni'
  read(*,*) e
  open(1,file='basis.dat')
  nso=e*2
  
  max2=2**0
  do i=2,13
     max2=max2+2**i
    
  enddo
    max=max2+2**25
  write(1,*)max
  min=0
  do i=0,13
     min=min+2**i
  enddo
  write(1,*)min
  write(2,*) n
endprogram basis
