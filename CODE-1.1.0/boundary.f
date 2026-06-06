      subroutine boundary(xv,il,iu,jl,ju,lls,nvp)
      use gridsetup
      use relax
      use xve
      use metryic
      use msga

      Implicit None

      integer,intent(in) :: il,iu,jl,ju,lls,nvp
      real xv(il:iu,jl:ju,lls,nvp)
      integer :: i,j,k,nvs !,ia,ja
      real :: speed1,speed2    

	  do 20 nvs=1,nv
	  do 20 k=1,l
      do 20 j=1,mp
      do 20 i=1,np
	                xv(i,j,k,nvs)=xv(i,j,k,nvs)*
     &	        (1.-relaxxv(i,j,k,min(nvs,4)))+xe(i,j,k,nvs)*relaxxv(i,j,k,min(nvs,4))
   20 continue

      if(ibctopbot.eq.0) then ! generally 1 for firetec
        do j=1,mp
        do i=1,np
         if(islip.eq.1) xv(i,j,1,1)=0.
         if(islip.eq.1) xv(i,j,1,2)=0.
         speed1=sqrt((xv(i,j,1,1)/xv(i,j,1,nv))**2
     &           +(xv(i,j,1,2)/xv(i,j,1,nv))**2)
         speed2=sqrt((xv(i,j,2,1)/xv(i,j,2,nv))**2
     &           +(xv(i,j,2,2)/xv(i,j,2,nv))**2+.0000001)
         xv(i,j,1,3)=-(xv(i,j,1,1)*c13(i,j)*gmul(1)+xv(i,j,1,2)*c23(i,j)
     1                                  *gmul(1))/gi(i,j,1)
         xv(i,j,l,3)=0.
         if (iturb.ge.1) then
           xv(i,j,1,5)=xv(i,j,2,5)*min(speed2,speed1)/speed2
           xv(i,j,1,5)=0.0
           xv(i,j,1,6)=xv(i,j,2,6)*min(speed2,speed1)/speed2
         endif
        enddo
        enddo
      endif

      return
      end
