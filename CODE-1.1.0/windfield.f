c234567****************************************************
      module windfield

      implicit none


      save

      contains
c234567****************************************************
        subroutine windfld_frqwriteio(fname)

      use gridsetup
      use xvo
      use pres
      use msga

      implicit none

      character(len=60):: fname

         if(mpi_rank.eq.0)then
           open(unit=21,file=fname,form='unformatted',
     +          status='unknown')
         endif
            call windfld_writeio(xvb(1-ih,1-ih,1,1),21,
     .1-ih,np+ih,1-ih,mp+ih,l)
            call windfld_writeio(xvb(1-ih,1-ih,1,2),21,
     .1-ih,np+ih,1-ih,mp+ih,l)
            call windfld_writeio(xvb(1-ih,1-ih,1,3),21,
     .1-ih,np+ih,1-ih,mp+ih,l)
            call windfld_writeio(xvb(1-ih,1-ih,1,4),21,
     .1-ih,np+ih,1-ih,mp+ih,l)
            if(iturb.ge.1)
     .      call windfld_writeio(xvb(1-ih,1-ih,1,5),21,
     .1-ih,np+ih,1-ih,mp+ih,l)
            if(iturb.ge.1)
     .      call windfld_writeio(xvb(1-ih,1-ih,1,6),21,
     .1-ih,np+ih,1-ih,mp+ih,l)
            call windfld_writeio(xvb(1-ih,1-ih,1,7),21,
     .1-ih,np+ih,1-ih,mp+ih,l)
            call windfld_writeio(xvb(1-ih,1-ih,1,nv),21,
     .1-ih,np+ih,1-ih,mp+ih,l)
c            call windfld_writeio(pr(1-ih,1-ih,1),21,
c     .1-ih,np+ih,1-ih,mp+ih,l)
         if(mpi_rank.eq.0)then
            close(21)
         endif
      return
        end subroutine windfld_frqwriteio

c***********************************************************
        subroutine windfld_writeio(data,iunit
     &,il,iu,jl,ju,nzdim)
      use gridsetup
      use msga

      Implicit None

c
      integer  il,iu,jl,ju,nzdim,iunit
      integer  i,j,k,iproc,iprocx,jprocy,ia,ja,iia,jja,nsize
      real data(il:iu,jl:ju,nzdim)
c
      real :: chtemp
c
      real,allocatable::tmparray(:,:,:,:),outdata(:,:,:)
     &,outdataj(:,:,:)
c
      if (mpi_rank.eq.0)
     +  allocate(tmparray(il:iu,jl:ju,nzdim,nproc),
c     +            outdata(ibcells*2,m,nzdim),
     +            outdata(ibcells*2,je-js+1,nzdim),
     +            outdataj(ie-is+1,jbcells*2,nzdim))
c     +            outdataj(n,jbcells*2,nzdim))
        nsize=(iu-il+1)*(ju-jl+1)*nzdim
        call mpi_gather(data,nsize,mpi_real,tmparray,nsize,mpi_real,
     +     0,mpi_comm_world,ierror)
      if (mpi_rank.eq.0) then
        do iprocx=1,nprocx
          do jprocy=1,nprocy
           iproc=1+(iprocx-1)+(jprocy-1)*nprocx
           do k=1,nzdim
             do j=1,mp
               do i=1,np
                 ia=(iprocx-1)*np + i
                 ja=(jprocy-1)*mp + j
                 if (ia.ge.is.and.ia.le.ie.and.ja.ge.js.and.ja.le.je) then
                 if(ia.le.(ibcells+is-1))then
                   iia=ia-is+1
                   jja=ja-js+1
                   outdata(iia,jja,k)=tmparray(i,j,k,iproc)
                   !write(6,*) 'out:iia,jja=',iia,jja
                 endif
                 if(ia.gt.ie-ibcells)then
                   iia=ia-(ie-ibcells*2)
                   jja=ja-js+1
                   outdata(iia,jja,k)=tmparray(i,j,k,iproc)
                   !write(6,*) 'out:iia,jja=',iia,jja
                 endif
                 if(ja.le.(jbcells+js-1))then
                   iia=ia-is+1
                   jja=ja-js+1
                   outdataj(iia,jja,k)=tmparray(i,j,k,iproc)
                   !write(6,*) 'outj:iia,jja=',iia,jja
                 endif
                 if(ja.gt.je-jbcells)then
                   iia=ia-is+1
                   jja=ja-(je-jbcells*2)
                   outdataj(iia,jja,k)=tmparray(i,j,k,iproc)
                   !write(6,*) 'outj:iia,jja=',iia,jja
                  endif
                 endif
               enddo
             enddo
           enddo
         enddo
        enddo
c
c       output file
        chtemp=sum(outdata)
        chtemp=sum(outdata)
        write (iunit)outdata,outdataj
c       write(6,*)
c    &'WINDFIELD WRITE ',outdata(ibcells,n/2,5),outdata(ibcells*2,n/2,5),
c    &outdataj(m/2,jbcells,5),outdataj(m/2,jbcells*2,5)
c
        deallocate(tmparray,outdata,outdataj)
      endif
c
      return

        end subroutine windfld_writeio

c234567***************************************
c        subroutine windfld_read(fname,xvbdata,prdata)
        subroutine windfld_read(fname,xvbdata)

      use gridsetup
      use xvbin
      use msga

      implicit none


      character(len=60):: fname
      real xvbdata(1-ih:np+ih,1-ih:mp+ih,l,nv)
      !real prdata(1-ih:np+ih,1-ih:mp+ih,l)

        if(mpi_rank.eq.0)then
          open(41,file=fname,form='unformatted',status='old')
        endif
        call    windfld_readio(xvbdata(1-ih,1-ih,1,1),41,l,!1,
     .1-ih,np+ih,1-ih,mp+ih)
        call    windfld_readio(xvbdata(1-ih,1-ih,1,2),41,l,!2,
     .1-ih,np+ih,1-ih,mp+ih)
        call    windfld_readio(xvbdata(1-ih,1-ih,1,3),41,l,!3,
     .1-ih,np+ih,1-ih,mp+ih)
        call    windfld_readio(xvbdata(1-ih,1-ih,1,4),41,l,!4,
     .1-ih,np+ih,1-ih,mp+ih)
        call    windfld_readio(xvbdata(1-ih,1-ih,1,5),41,l,!5,
     .1-ih,np+ih,1-ih,mp+ih)
        call    windfld_readio(xvbdata(1-ih,1-ih,1,6),41,l,!6,
     .1-ih,np+ih,1-ih,mp+ih)
        call    windfld_readio(xvbdata(1-ih,1-ih,1,7),41,l,!7,
     .1-ih,np+ih,1-ih,mp+ih)
        call    windfld_readio(xvbdata(1-ih,1-ih,1,8),41,l,!8,
     .1-ih,np+ih,1-ih,mp+ih)
c        call    windfld_readio(prdata(1-ih,1-ih,1),41,l,9,
c     .1-ih,np+ih,1-ih,mp+ih)
        if(mpi_rank.eq.0)then
          close(41)
        endif

        return
        end subroutine windfld_read

c***********************************************************
        subroutine windfld_readio(array,iunit,nzdim,il,iu,jl,ju) ! (array,iunit,nzdim,id,il,iu,jl,ju)
      use gridsetup
      use msga
c
      Implicit None
c

      integer  il,iu,jl,ju,nzdim,iunit
      integer  i,j,k,iproc,iprocx,jprocy,ia,ja,ii,jj,nsize !id
      real array(np,mp,nzdim)
      real,allocatable:: tmparray(:,:,:,:),indata(:,:,:),
     &indataj(:,:,:),indata2(:,:,:)
c
      real :: indatamax,indatamin
c

      if (mpi_rank.eq.0) then
        allocate(indata(ibcells*2,m,nzdim),indata2(n,m,nzdim),
     &           indataj(n,jbcells*2,nzdim),
     +           tmparray(il:iu,jl:ju,nzdim,nproc))
c    +           tmparray(np,mp,nzdim,nproc))
        read (iunit)indata,indataj

        indatamax=0.0
        indatamin=0.0

        do k=1,nzdim
          do j=1,m
            do i=1,n
                 if(i.le.ibcells)then
                   indata2(i,j,k)=indata(i,j,k)
                 elseif(i.gt.n-ibcells)then
                   ii=i-(n-ibcells*2)
                   indata2(i,j,k)=indata(ii,j,k)
                 else
                   indata2(i,j,k)=0.0
                 endif
                 if(j.le.jbcells)then
                   indata2(i,j,k)=indataj(i,j,k)
                 elseif(j.gt.m-jbcells)then
                   jj=j-(m-jbcells*2)
                   indata2(i,j,k)=indataj(i,jj,k)
                 endif
                 indatamax=max(indata2(i,j,k),indatamax)
                 indatamin=min(indata2(i,j,k),indatamin)
            enddo
          enddo
        enddo

c       chtmp=sum(indata)

        if(mpi_rank.eq.0)then
c       print *,' JMC READIO ',id,indata(5,5,5),indatamax,indatamin
        endif
        do iprocx=1,nprocx
          do jprocy=1,nprocy
           iproc=1+(iprocx-1)+(jprocy-1)*nprocx
           do k=1,nzdim
             do j=1,mp
               do i=1,np
                 ia=(iprocx-1)*np + i
                 ja=(jprocy-1)*mp + j
                 if(ia.le.ibcells)then
                   tmparray(i,j,k,iproc)=indata2(ia,ja,k)
                 elseif(ia.gt.n-ibcells)then
                   tmparray(i,j,k,iproc)=indata2(ia,ja,k)
                 else
                   tmparray(i,j,k,iproc)=0.0
                 endif
                 if(ja.le.jbcells)then
                   tmparray(i,j,k,iproc)=indata2(ia,ja,k)
                 elseif(ja.gt.m-jbcells)then
                   tmparray(i,j,k,iproc)=indata2(ia,ja,k)
                 endif
               enddo
             enddo
           enddo
         enddo
        enddo
      endif

      nsize=(iu-il+1)*(ju-jl+1)*nzdim
c     nsize=(np)*(mp)*nzdim
      call mpi_scatter(tmparray,nsize,mpi_real,array,nsize,mpi_real,
     +   0,mpi_comm_world,ierror)
      if(mpi_rank.eq.0)deallocate(tmparray,indata,indata2)

      return
       end subroutine windfld_readio

      end module windfield
