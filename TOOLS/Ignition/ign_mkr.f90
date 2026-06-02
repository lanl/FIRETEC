!!!!!!!*******************************************************************!!!!!
!      ign_mkr.f90: The purpose of this piece of code will be to generate an
!      ignition file which will be read into firetec to ignite the fire in
!      a manner specified by the parameters contained within the file.
!      Author: Jeremy Sauer <jsauer@lanl.gov>  
!      Date: 1/13/06
!      inputs:
!      outputs: ignite.dat
!!!!!!!*******************************************************************!!!!!
 
      Program ign_mkr
      Implicit None
     
      integer::i,j,k,icnt 
      integer::igntype,numignpts,numfires
      integer,allocatable::ignpts(:,:),ila(:),iua(:),jla(:),jua(:)
      real::targettemp,startigntime,endigntime,ramprate
      real::xfirelinelow,xfirelinehigh,yfirelinelow,yfirelinehigh
      real::flamedistance
      integer:: il,iu,jl,ju,numi,numj,tmp

      open(unit=1000,file='ignite.dat',&
      form='formatted',status='unknown')
 
      write(6,*) "Please input an ignition type & 
      (0=none,1=rinitfire,2=terratorch,3=multiple fires):" 
      read(5,'(I1)') igntype
 
      if(igntype==1)then     !write out a rinitfire style ignition file!

       write(6,*) "What is the target temperature (K) & 
       you would like to reach?(i.e. 1000.0)"
       targettemp = 1000.0
       read(5,*) targettemp
       write(6,*) "What time (s) would you like &
       to start ignition?(i.e. 0.0)"
       startigntime=0.0
       read(5,*) startigntime 
       write(6,*) "At what rate (K/s) would you &
       like to ramp up the ignition &
       cell temperature?(i.e. 350.0)"
       ramprate = 350.0
       read(5,*) ramprate 
       write(6,*) "Please specify a lower &
       bound cell index in the x-direction."
       read(5,'(I5)') il
       write(6,*) "Please specify a upper bound &
       cell index in the x-direction."
       read(5,'(I5)') iu
       write(6,*) "Please specify a lower &
      bound cell index in the y-direction."
       read(5,'(I5)') jl
       write(6,*) "Please specify a upper bound &
       cell index in the y-direction."
       read(5,'(I5)') ju
     
       numi = iu-il+1
       numj = ju-jl+1
       numignpts = numi*numj 
       allocate(ignpts(numignpts,2))


       i = il
       j = jl
       do tmp=1,numignpts
        ignpts(tmp,1)=i
        ignpts(tmp,2)=j
        if(j==ju)then
         i=i+1
         j=jl
        else
         j=j+1
        endif
       enddo
        
       write(1000,'(A15,I5)') 'igntype=',igntype
       write(1000,'(A15)') '&rinitlist'       
       write(1000,'(A15,I5)') 'numignpts=',numignpts
       write(1000,'(A15,F8.2)') 'targettemp=',targettemp
       write(1000,'(A15,F8.2)') 'startigntime=',startigntime
       write(1000,'(A15,F8.2)') 'ramprate=',ramprate
       write(1000,'(A1)') '/'
       do i=1,numignpts
        write(1000,'(2I5)') ignpts(i,1),ignpts(i,2) 
       enddo

       deallocate(ignpts)

      elseif(igntype==2)then  !write out a terratorch style ignition file!

       write(6,*) "What is the target temperature (K) &
       you would like to reach?"
       read(5,*) targettemp
       write(6,*) "What time (s) would you like to start ignition?"
       read(5,*) startigntime
       write(6,*) "What time (s) would you like to end ignition?"
       read(5,*) endigntime
       write(6,*) "What is the lower bound coordinate (m) &
       in the x direction?"
       read(5,*) xfirelinelow 
       write(6,*) "What is the upper bound coordinate (m) &
       in the x direction?"
       read(5,*) xfirelinehigh 
       write(6,*) "What is the lower bound coordinate (m) &
       in the y direction?"
       read(5,*) yfirelinelow 
       write(6,*) "What is the upper bound coordinate (m) &
       in the y direction?"
       read(5,*) yfirelinehigh 
       write(6,*) "What is flame halo distance (m)?"
       read(5,*) flamedistance 

       write(1000,'(A15,I5)') 'igntype=',igntype
       write(1000,'(A15)') '&torchlist'      
       write(1000,'(A15,F8.2)') 'targettemp=',targettemp
       write(1000,'(A15,F8.2)') 'startigntime=',startigntime
       write(1000,'(A15,F8.2)') 'endigntime=',endigntime
       write(1000,'(A15,F8.2)') 'xfirelinelow=',xfirelinelow
       write(1000,'(A15,F8.2)') 'xfirelinehigh=',xfirelinehigh
       write(1000,'(A15,F8.2)') 'yfirelinelow=',yfirelinelow
       write(1000,'(A15,F8.2)') 'yfirelinehigh=',yfirelinehigh
       write(1000,'(A15,F8.2)') 'flamedistance=',flamedistance
       write(1000,'(A1)') '/'

      elseif(igntype==3)then     !JMC write out a multiple fire style ignition file!

       write(6,*) "What is the target temperature (K) &
      you would like to&
      reach?(i.e. 1000.0)"
       targettemp = 1000.0
       read(5,*) targettemp
       write(6,*) "What time (s) would you &
       like to start ignition?(i.e. 0.0)"
       startigntime=0.0
       read(5,*) startigntime 
       write(6,*) "At what rate (K/s) would you like to ramp up the  &
       ignition cell temperature?(i.e. 350.0)"
       ramprate = 350.0
       read(5,*) ramprate 
       write(6,*) "How many fires would you like to start?(i.e. 2)"
       numfires=2
       read(5,*) numfires
       allocate(ila(numfires),iua(numfires),jla(numfires),jua(numfires))
       numignpts=0
       do k=1,numfires !JMC loop through the number of fires to start
       write(6,*) "Please specify a lower bound cell index in &
        x-direction",k
       read(5,'(I5)') ila(k)
       write(6,*) "Please specify a upper bound cell index in &
       x-direction",k
       read(5,'(I5)') iua(k)
       write(6,*) "Please specify a lower bound cell index in &
       y-direction",k
       read(5,'(I5)') jla(k)
       write(6,*) "Please specify a upper bound cell index in &
       y-direction",k
       read(5,'(I5)') jua(k)
     
       numi = iua(k)-ila(k)+1
       numj = jua(k)-jla(k)+1
       numignpts = numi*numj + numignpts
       enddo
       allocate(ignpts(numignpts,2))


       icnt=1
       do
        do k=1,numfires
         do i=ila(k),iua(k)
          do j=jla(k),jua(k)
           ignpts(icnt,1)=i
           ignpts(icnt,2)=j
           icnt=icnt+1
          enddo
         enddo
        enddo
        if(icnt>numignpts)exit
       enddo
        
       write(1000,'(A15,I5)') 'igntype=',igntype
       write(1000,'(A15)') '&rinitlist'       
       write(1000,'(A15,I5)') 'numignpts=',numignpts
       write(1000,'(A15,F8.2)') 'targettemp=',targettemp
       write(1000,'(A15,F8.2)') 'startigntime=',startigntime
       write(1000,'(A15,F8.2)') 'ramprate=',ramprate
       write(1000,'(A1)') '/'
       do i=1,numignpts
        write(1000,'(2I5)') ignpts(i,1),ignpts(i,2) 
       enddo

       deallocate(ila,iua,jla,jua)
       deallocate(ignpts)

      elseif(igntype==0)then  !no fire!
       write(1000,'(A15,I5)') 'igntype=',0
      endif

 
      close(1000) 
    
      write(6,*) "Thank for your inputs!" 
      write(6,*) "Your new ignition file is ready as ignite.dat!"    
      End Program ign_mkr
