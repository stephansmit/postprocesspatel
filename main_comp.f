c     spectral in theta and z, finite difference in r.
c
c
c     15 - 12- 2009    full skew symmetric terms including divergence. 
c                        Non uniform grid in r (routine mkgrid)
c                        Adamsb-Bashforth 3 time integration 
c
c
c     28-2-2012   parallel version with 2D parallization using 2decomp http://www.2decomp.org
c
c
c     Bendiks Jan Boersma, Delft University of Technology, Delft, The Netherlands.
c
c
     

      use decomp_2d
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'

      integer ini,iload,ploc,nstap,ierr,istap,ini_sol,icount,i_dex,j_dex
      real ut(jmax),dut(jmax),ener,t_som,bulk_tot, stress_tot,f_colebrook,f_shokling
      real Wbulk, ReB,bulk,str1,str2,time,dpdx,Lz,Lt,ti_s(imax),ti_c(0:imax),vm(imax)
      real cbulk,Twh,Twc,Qw,fl1,fl2
      integer cbc

      call mpi_init(ierr)
      call decomp_2d_init(imax+2,jmax,kmax,p_row,p_col)
      rank = nrank
      if (nrank.eq.0) then
      write(6,*) xsize(1),xsize(2),xsize(3),xsize(1)*xsize(2)*xsize(3)
      write(6,*) ysize(1),ysize(2),ysize(3),ysize(1)*ysize(2)*ysize(3)
      write(6,*) zsize(1),zsize(2),zsize(3),zsize(1)*zsize(2)*zsize(3)
      endif

      Re = 395.
      Pr = 1.
      Lz = 8.*atan(1.)
      Lt = 4.*atan(1.)
      Twh=1.
      Twc=1. 
      cbc=0 !0 for isotemp... 1 for isoflux
      Twb=Twh
      Twt=Twc

      dt =1e-4
      nstap=300000

      dpdx= 1.0


      iload = 1
  



      open(11,file='icount')
      read(11,*)icount,time
      close(11)
   
      call mkgrid(Lz,Lt,nrank)
      if (iload.eq. 0)   call init (nrank)
      if (iload.ne. 0)   call loadd(0,nrank,icount)

      visc=1./Re
      lambdanew=cnew/(Re*Pr)
      time = 0.
      ini = 0
      ini_sol=0

      call bound(unew,vnew,wnew,nrank)
      call bounds(cnew,Twh,Twc,qw,nrank)
      call output(0,rank)
      call chkdt(nrank)
      if (nrank.eq.0) call report
      do istap=1,nstap
      if (mod(istap,100).eq. 0) then
        call cmpbs(bulk,cbulk,str1,str2,fl1,fl2)
  
      if (nrank.eq.0) write(6,111) time, RE*bulk , str1, str2, dpdx
111   format ('tijd= ', f16.5, ' RE_bulk= ',F16.6, ' Stress1= ',F16.6, ' Stress2= ',F16.6, ' dp/dx= ', F16.6)
      if (nrank.eq.0) write(6,200) dt, cbulk, fl1, fl2, istap
200   format (' dt= ', f16.6, ' cbulk= ', f16.6, ' flux1= ',F16.6, ' flux2= ', F16.6, ' istap=', I6.6)
      stime = MPI_WTIME()
      endif

      call adamsb(ini,dpdx,nrank)
      lambdanew=cnew/(Re*Pr)
      ini = 1

      if (mod(istap,100).eq. 0) then
      if (rank.eq.0) write(6,*) ' Adamsb tijd = ', MPI_WTIME()-stime
      stime = MPI_WTIME()
      endif

      call bound(dudt,dvdt,dwdt,rank)
      call bounds(dcdt,Twh,Twc,Qw,rank)
      call trunc(dudt,rank)
      call trunc(dvdt,rank)
      call trunc(dwdt,rank)
      call trunc(dcdt,rank)

      if (mod(istap,100).eq. 0) then
      if (rank.eq.0) write(6,*) ' trunc  tijd = ', MPI_WTIME()-stime
      stime = MPI_WTIME()
      endif

      call fillps

      if (mod(istap,100).eq. 0) then
      if (rank.eq.0) write(6,*) ' fillps tijd = ', MPI_WTIME()-stime
      stime = MPI_WTIME()
      endif
     
      if (ipois.eq.2) then
         call solver_2(ini_sol,p,ru,rp,dr,dtheta,dz,mr_s,mr_c,nrank)
      endif
      if (ipois.ne.2) then 
      call solver_com(ini_sol,p,ru,rp,dr,dtheta,dz,mr_s,mr_c,nrank)
      ini_sol =10
      endif

      if (mod(istap,100).eq. 0) then
      if (rank.eq.0) write(6,*) ' Solver tijd = ', MPI_WTIME()-stime
      stime = MPI_WTIME()
      endif

      call correc

      if (mod(istap,100).eq. 0) then
      if (rank.eq.0) write(6,*) ' correc tijd = ', MPI_WTIME()-stime
      endif

      call bound(unew,vnew,wnew,rank)
      call bounds(cnew,Twh,Twc,Qw,rank)

      if (mod(istap,100).eq. 0) call chkdiv(rank)
      if (mod(istap,100).eq. 0) call chkdt (rank)

      if (rank.eq.0     .and.  mod(istap,1000).eq.0  ) write(6,*) ' tijd     =    ',time,dt
      if (rank.eq.0     .and.  mod(istap, 200).eq. 0) write(6,456) ener(),stress,bulk,dpdx
      if (rank.eq.0     .and.  mod(istap, 200).eq. 0) write(16,456) ener(),stress,bulk,dpdx
      if (rank.eq.0     .and.  mod(istap, 200).eq. 0) call Etotexp(time)
456   format( ' The turbulent energy | stress | bulk | dp/dx  = ',F13.6 ,' | ', F13.6, ' | ', F13.6, '|', F13.6)
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if (mod(istap,1000).eq.0) call output(istap/500,rank)
      if (mod(istap,1000).eq.0) then
      icount = icount + 1
      call loadd(1,nrank,icount)
      if (nrank.eq.0) then
      open(11,file='icount')
      write(11,*) icount,time
      close(11)
      endif
      endif
      
      if (istap.eq. 999999) call output2(istap/1500,rank)

      time = time + dt 
      enddo
      call decomp_2d_finalize
      call mpi_finalize(ierr)
      stop
      end
 

      subroutine loadd(ini,nnrank,istap)
      use decomp_2d
      use decomp_2d_io
      include 'param.txt'
      include 'common.txt'
      real utmp(0:i1,jmax/p_row,kmax/p_col)
      real vtmp(0:i1,jmax/p_row,kmax/p_col)
      real wtmp(0:i1,jmax/p_row,kmax/p_col)
      real ptmp(0:i1,jmax/p_row,kmax/p_col)
      real dum1(imax),dum2(imax)

      integer istap
      character*5 cha
      character*5 cha2
      call cnvstr(nnrank,cha)
      call cnvstr(istap,cha2)

      do i=1,imax
       ptmp(i,:,:)=p(i,:,:)
      enddo

      if (ini.eq.0) then
      call decomp_2d_read_one(1,unew,'DNS/u.'//cha2//'.dns')
      call decomp_2d_read_one(1,vnew,'DNS/v.'//cha2//'.dns')
      call decomp_2d_read_one(1,wnew,'DNS/w.'//cha2//'.dns')
      call decomp_2d_read_one(1,cnew,'DNS/c.'//cha2//'.dns')
!      cnew=0.995*cnew

!     open(19,file = 'temp')
!     do i=1,imax
!       read(19,*) dum1(i),dum2(i)
!     enddo
!     close(19)
!     
!     do i=1,imax
!       cnew(i,:,:)=dum2(i)+(cnew(i,:,:)-dum1(i))/sqrt(cnew(i,:,:))
!     enddo

!     open(19,file = 'mean')
!     do i=1,imax
!       read(19,*) dum1(i),dum2(i)
!     enddo
!     close(19)
!
!     do i=1,imax
!       wnew(i,:,:)=dum2(i)+(wnew(i,:,:)-dum1(i))/sqrt(cnew(i,:,:))
!       unew(i,:,:)=unew(i,:,:)/sqrt(cnew(i,:,:))
!       vnew(i,:,:)=vnew(i,:,:)/sqrt(cnew(i,:,:))
!     enddo



      endif

      if (ini.eq.1) then
      call decomp_2d_write_one(1,unew,'DNS/u.'//cha2//'.dns')
      call decomp_2d_write_one(1,vnew,'DNS/v.'//cha2//'.dns')
      call decomp_2d_write_one(1,wnew,'DNS/w.'//cha2//'.dns')
      call decomp_2d_write_one(1,cnew,'DNS/c.'//cha2//'.dns')
      call decomp_2d_write_one(1,ptmp,'DNS/p.'//cha2//'.dns')
      endif


      end



      subroutine output(istap,rank)
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer istap
      character*5 cha
      character*5 cha2
      real um(imax),vm(imax),wm(imax),cm(imax)
      real umm(imax),vmm(imax),wmm(imax),cmm(imax),tcmm(imax)
      real ur(imax),vr(imax),wr(imax),cr(imax),hi_s(imax),ci_s(imax),hi_c(0:imax),ci_c(0:imax)
      real uv(imax),uw(imax),vw(imax),cu(imax),cv(imax),cw(imax),tcm(imax)
      um=0
      vm=0
      wm=0
      cm=0
      ur=0
      vr=0
      wr=0
      cr=0
      uw=0
      uv=0
      vw=0
      cu=0
      cv=0
      cw=0
      do k=1,kmax/p_col
	 do j=1,jmax/p_row
	    do i=1,imax
	      um(i)=um(i)+unew(i,j,k)/(jmax*kmax)
	      vm(i)=vm(i)+vnew(i,j,k)/(jmax*kmax)
	      wm(i)=wm(i)+wnew(i,j,k)/(jmax*kmax)
              cm(i)=cm(i)+cnew(i,j,k)/(jmax*kmax)
             tcm(i)=tcm(i)+lambdanew(i,j,k)/(jmax*kmax)
	    enddo
	 enddo
      enddo
      call mpi_allreduce(um,umm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(vm,vmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(wm,wmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(cm,cmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(tcm,tcmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier) 
      um = umm
      vm = vmm
      wm = wmm
      cm = cmm
      tcm=tcmm
      do k=1,kmax/p_col
	do j=1,jmax/p_row
	  do i=1,imax
	    ur(i)=ur(i)+(unew(i,j,k)-um(i))**2
	    vr(i)=vr(i)+(vnew(i,j,k)-vm(i))**2
	    wr(i)=wr(i)+(wnew(i,j,k)-wm(i))**2
            cr(i)=cr(i)+(cnew(i,j,k)-cm(i))**2
            vw(i)=vw(i)+(vnew(i,j,k)-vm(i))*(wnew(i,j,k)-wm(i))
            cv(i)=cv(i)+(vnew(i,j,k)-vm(i))*(cnew(i,j,k)-cm(i))
            cw(i)=cw(i)+(cnew(i,j,k)-cm(i))*(wnew(i,j,k)-wm(i))
	  enddo
          do i=1,imax
            uv(i)=uv(i)+(unew(i,j,k)+unew(i-1,j,k)-um(i)-um(i-1))*(vnew(i,j,k)-vm(i))*0.5
            uw(i)=uw(i)+(unew(i,j,k)+unew(i-1,j,k)-um(i)-um(i-1))*(wnew(i,j,k)-wm(i))*0.5
            cu(i)=cu(i)+(unew(i,j,k)+unew(i-1,j,k)-um(i)-um(i-1))*(cnew(i,j,k)-cm(i))*0.5
          enddo
	enddo
      enddo
      call mpi_allreduce(ur,umm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(vr,vmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(wr,wmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(cr,cmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ur = umm
      vr = vmm
      wr = wmm
      cr = cmm
      call mpi_allreduce(uw,wmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(cu,cmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      uw = wmm 
      cu = cmm
      hi_s = wm
      call der1w_s_c6(imax,hi_s,hi_c,dr)
      hi_c = hi_c * mr_c
      hi_s = cm
      call der1c_s_c6(imax,hi_s,ci_c,dr,Twb,Twt)
      ci_c = ci_c * mr_c
      call inter_c_s(imax,ci_c,ci_s,dr)
      call inter_c_s(imax,hi_c,hi_s,dr)
      if (rank.eq.0) then
      open(17,file = 'prof')
      open(18,file = 'prof_log')
      open(19,file = 'prof_temp')
      do i=1,imax
	write(17,112) rp(i),um(i),vm(i),wm(i),
     1 sqrt(ur(i)/(jmax*kmax)),sqrt(vr(i)/(jmax*kmax)),sqrt(wr(i)/(jmax*kmax)),(uw(i)/(jmax*kmax)),
     2 -visc*hi_s(i),(uw(i)/(jmax*kmax))-visc*hi_s(i)
      enddo
      do i=1,imax
        write(19,112) rp(i),cm(i),sqrt(cr(i)/(jmax*kmax)),(cu(i)/(jmax*kmax)),
     2 -tcm(i)*ci_s(i),(cu(i)/(jmax*kmax))-tcm(i)*ci_s(i),tcm(i)
      enddo
      close(17)
      close(18)
      close(19)
112   FORMAT(10E16.5) 
      endif
      end



      subroutine fillps
c
c  right hand side of the poisson equation
c
c       1 dru*    1 dv*    dw* 
c       - --- +   - ---  + --- = 0
c       r d r     r d t    d z
c
c
c         (1)      (2)      (3)
c
      use decomp_2d
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer im,ier,idex,i_dex
      real ft(jmax),dft(jmax),fz(kmax),dfz(kmax),som,stime
      real hi_s(imax,jmax/p_row,kmax/p_col),hi_c(0:imax,jmax/p_row,kmax/p_col)
      real tmp (0:imx,jmax,kmax/p_col)
      real tmp2(0:imx,jmax,kmax/p_col)
      real tmp3(0:imx,jmax/p_col,kmax)
      real p1(0:imax+1,jmax/p_row,kmax/p_col)
      p = 0
      

      call transpose_x_to_y(dvdt,tmp) 

      do i=0,imx
	do k=1,kmax/p_col
	  do j=1,jmax
	    ft(j)=tmp(i,j,k)
	  enddo
	  call four1(jmax,ft,dft,dtheta)
	  do j=1,jmax
	   tmp(i,j,k)=dft(j)
	  enddo
	enddo
      enddo
      call transpose_x_to_y(dwdt,tmp2)
      call transpose_y_to_z(tmp2,tmp3)
c term (1)
      do i=0,imx
	 do j=1,jmax/p_col
	   do k=1,kmax
	     fz(k)=tmp3(i,j,k)
	   enddo
	   call four1(kmax,fz,dfz,dz)
	   do k=1,kmax
	    tmp3(i,j,k)=dfz(k)
	   enddo
	 enddo
	enddo
        call transpose_z_to_y(tmp3,tmp2)
        tmp = tmp + tmp2 
        call transpose_y_to_x(tmp,p1)
       
        if (ipois.ne.2) then 
         do k=1,kmax/p_col
	  do j=1,jmax/p_row
	   do i=0,imax
            hi_c(i,j,k)=dudt(i,j,k)
           enddo
          enddo
         enddo
           call deriv_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr)
         do k=1,kmax/p_col
	    do j=1,jmax/p_row
           do i=1,imax
	     p(i,j,k)=p1(i,j,k)+hi_s(i,j,k)*mr_s(i)
           enddo
         enddo
         enddo
       endif

       if (ipois.eq.2) then

         do k=1,kmax/p_col
            do j=1,jmax/p_row
           do i=1,imax
         p(i,j,k)=p1(i,j,k)+(dudt(i,j,k)-dudt(i-1,j,k))/((ru(i)-ru(i-1)))
          enddo
         enddo
         enddo
      endif


       p = p /dt
       end
     
     
      subroutine chkdiv(rank)
      use decomp_2d
c    calculates divergence
c
c
c       1 dru     1 d v    d w 
c       - --- +   - ---  + --- = 0
c       r  dr     r dt     dz
c
c term    (1)      (2)     (3)
c

      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'

      integer im,ier,idex,i_dex,ierr
      real ft(jmax),dft(jmax),fz(kmax),dfz(kmax),som,dmaxx,div(0:i1,jmax/p_row,kmax/p_col)
      real t_dmaxx
      real t_som,hi_c(0:imax),hi_s(imax)
      real tmpx(0:i1 ,jmax/p_row,kmax/p_col)         
      real tmpy(0:imx,jmax      ,kmax/p_col)         
      real tmpz(0:imx,jmax/p_col,kmax      )         
      idex = i_dex(nrank)
      div = 0.
      dmaxx =0.
      call transpose_x_to_y(wnew,tmpy)
      call transpose_y_to_z(tmpy,tmpz)
      do i=0,imx
	 do j=1,jmax/p_col
	   do k=1,kmax
	    fz(k) =tmpz(i,j,k)
	   enddo
	   call four1(kmax,fz,dfz,dz)
	   do k=1,kmax
	    tmpz(i,j,k)=dfz(k)
	   enddo
	 enddo
      enddo
      call transpose_z_to_y(tmpz,tmpy)
      call transpose_y_to_x(tmpy,div)

      
      call transpose_x_to_y(vnew,tmpy)
      do i=0,imx
	do k=1,kmax/p_col
	  do j=1,jmax
	    ft(j)=tmpy(i,j,k)
	  enddo
	  call four1(jmax,ft,dft,dtheta)
	  do j=1,jmax
	   tmpy(i,j,k)=dft(j)
	  enddo
	enddo
      enddo
      call transpose_y_to_x(tmpy,tmpx)
      div = div  + tmpx
       som = 0 
c term (1)
       dmaxx = 0.
	   do j=1,jmax/p_row
	     do k=1,kmax/p_col
	 do i=0,imax
            hi_c(i)=unew(i,j,k)
         enddo
         call deriv_c_s(imax,hi_c,hi_s,dr)
         hi_s = hi_s*mr_s
          do i=1,imax
            if (ipois.ne.2)  div(i,j,k)=div(i,j,k)+hi_s(i)
            if (ipois.eq.2) div(i,j,k)=div(i,j,k)+(unew(i,j,k)-unew(i-1,j,k))/((ru(i)-ru(i-1)))
	  enddo
	 enddo
       enddo
       dmaxx = 0
       som =0
       do k=1,kmax/p_col
         do j=1,jmax/p_row
	   do i=1,imax
             div(i,j,k)=div(i,j,k)*dr*dtheta*dz
             dmaxx = max (dmaxx,div(i,j,k))
             som = som + div(i,j,k)
           enddo
         enddo
       enddo
       call mpi_allreduce(som,t_som    ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
       call mpi_allreduce(dmaxx,t_dmaxx,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
       if (nrank.eq.0) write(6,*) 'divtot    =  ',t_som, 'divmax   = ',t_dmaxx
       if (abs(t_dmaxx).gt.1e-3) then
       call mpi_finalize(ierr)
       stop 'chkdiv'
       endif 
       end
     
     



      subroutine bound(u,v,w,rank)
      include 'param.txt'
      real u(0:i1,jmax/p_row,kmax/p_col)
      real v(0:i1,jmax/p_row,kmax/p_col)
      real w(0:i1,jmax/p_row,kmax/p_col)
      do k=1,kmax/p_col
	do j=1,jmax/p_row
	  u(imax,j,k)=0.
	  v(i1  ,j,k)=-v(imax,j,k)	
	  w(i1  ,j,k)=-w(imax,j,k)	
	  u(0   ,j,k)=0.
	  v(0   ,j,k)=-v(1,j,k)
	  w(0   ,j,k)=-w(1,j,k)
	enddo
       enddo
       end

      subroutine bounds(c,Twh,Twc,qw,rank)
      include 'param.txt'
      real c(0:i1,jmax/p_row,kmax/p_col)
      real Twh,Twc,qw
      do k=1,kmax/p_col
        do j=1,jmax/p_row
          c(i1  ,j,k)= 2*Twc-c(imax,j,k)      
          c(0   ,j,k)= 2*Twh-c(1,j,k)
        enddo
       enddo
       end



      subroutine correc
      use decomp_2d
      implicit none
      include 'param.txt'
      include 'common.txt'
      real som
      integer idex,i_dex
      real fz(kmax),dfz(kmax)
      real ft(jmax),dft(jmax)
      real hi_s(  imax,jmax/p_row,kmax/p_col)
      real hi_c(0:imax,jmax/p_row,kmax/p_col)

      real p1x(0:imax+1,jmax/p_row,kmax/p_col)
      real p1y(0:imx   ,jmax      ,kmax/p_col)
      real p1z(0:imx   ,jmax/p_col,kmax      )


      
      do k=1,kmax/p_col
	do j=1,jmax/p_row
	  do i=1,imax
            hi_s(i,j,k)=p(i,j,k)
          enddo
        enddo
      enddo
      call deriv_s_c_m(imax,jmax*kmax/px,hi_s,hi_c,dr)
      do k=1,kmax/p_col
        do j=1,jmax/p_row
          do i=1,imax
            if (ipois.ne.2) unew(i,j,k)=dudt(i,j,k)-hi_c(i,j,k)*dt*mr_c(i)
            if (ipois.eq.2) unew(i,j,k)=dudt(i,j,k)-dt*(p(i+1,j,k)-p(i,j,k))/(rp(i+1)-rp(i))
	  enddo
	enddo
      enddo 

      do k=1,kmax/p_col
	do j=1,jmax/p_row
	  do i=1,imax
	   p1x(i,j,k)=p(i,j,k)
          enddo
        enddo
      enddo

      call transpose_x_to_y(p1x,p1y)
      call transpose_y_to_z(p1y,p1z)
      do k=1,kmax/p_col
	do i=0,imx
	   do j=1,jmax
	     ft(j)=p1y(i,j,k)
	   enddo
	   call four1(jmax,ft,dft,dtheta)
           do j=1,jmax
	     p1y(i,j,k)=dft(j)
           enddo
        enddo
      enddo
      call transpose_y_to_x(p1y,p1x)

      do k=1,kmax/p_col
	 do j=1,jmax/p_row
           do i=1,imax 
	   vnew(i,j,k)=dvdt(i,j,k)-p1x(i,j,k)*dt
	   enddo
	enddo
      enddo 
      do k=1,kmax/p_col
         do j=1,jmax/p_row
	   do i=1,imax
	    p1x(i,j,k)=p(i,j,k)
           enddo
         enddo
      enddo

      do j=1,jmax/p_col
	do i=0,imx
	  do k=1,kmax
	    fz(k)=p1z(i,j,k)
	  enddo
	  call four1(kmax,fz,dfz,dz)
	  do k=1,kmax
	   p1z(i,j,k)=dfz(k)
	  enddo
	enddo
      enddo
      call transpose_z_to_y(p1z,p1y)
      call transpose_y_to_x(p1y,p1x)

      do k=1,kmax/p_col
	 do j=1,jmax/p_row
          do i=1,imax
	  wnew(i,j,k)=dwdt(i,j,k)-p1x(i,j,k)*dt
	  enddo
	enddo
      enddo
 
      cnew=dcdt

      end
      
      subroutine init(rank)
      include 'param.txt'
      include 'common.txt'
      real yplus
 
      unew= 0
      vnew= 0
      wnew= 15
      do k=1,kmax/p_col
        do j=1,jmax/p_row
          do i=1,imax
!              yplus=(1.0-rp(i))/visc
!             if  (yplus .lt. 11.6) wnew(i,j,k)=yplus
!             if  (yplus .gt. 11.6)  wnew(i,j,k)= 2.5*log(yplus)+5.5
!            wnew(i,j,k)=wnew(i,j,k)+cos(1.*(rank+1)*i*j*k)*4
!            vnew(i,j,k)=sin(1.*j*j*k) 
            wnew(i,j,k)=(Re/2.)*(1.- (rp(i)-1.)**2.)
            cnew(i,j,k)=3./2.*((rp(i)-1.)**2./12.-(rp(i)-1.)**2./2.)+13./8.
            enddo
          enddo
       enddo

      end  

      subroutine mkgrid(Lz,Lt,rank)
      include 'param.txt'
      include 'common.txt'
      real const,Lz,Lt,rp_c(0:imax),rp_s(imax),delta(imax)
      dz =Lz/(kmax)
      dr = 2.0/imax
      dtheta =Lt/jmax
       rmax = 2.0
       ru(0)=0. 

       do i=1,imax/2
          x  = 1.*i/imax
          dx = 0.5-1.51*(x-0.5)**2.!0.00145*x**3. -5.966*x**2. +734.37*x+1493.81
          ru(i)=ru(i-1)+dx
        enddo
        rnorm = ru(imax/2)
        do i=1,imax/2
          ru(i)=ru(i)/rnorm
        enddo
       do i=imax,imax/2+1,-1
         ru(i)=2.-ru(imax-i)
       enddo


      do i=1,imax
        rp(i)=0.5*(ru(i)+ru(i-1))
        delta(i)=ru(i)-ru(i-1)
      enddo

      rp(0 )=ru(0)-rp(1)
      rp(i1)=ru(imax)+(Ru(imax)-rp(imax))
      do i=1,imax
	  rp_s(i)=rp(i)
      enddo
      call deriv_s_c(imax,rp_s,mr_c,dr)
      call inter_c_s(imax,mr_c,mr_s,dr)
      mr_c = 1./mr_c
      mr_s = 1./mr_s
      if (rank.eq.0) then
      open(11,file = 'grid.txt')
      write(11,*) Re,Ru(imax)
      do i=1,imax
         write(11,'(i5,4F12.6)') i,Ru(i),Rp(i),delta(i),delta(i)*Re
      enddo
      endif
      close(11)
      end 


      subroutine adamsb(ini,dpdx,rank)
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real dfr_n(0:imax+1,jmax/p_row,kmax/p_col)
      real dft_n(0:imax+1,jmax/p_row,kmax/p_col)
      real dfz_n(0:imax+1,jmax/p_row,kmax/p_col)
      real dfc_n(0:imax+1,jmax/p_row,kmax/p_col)
      real dfr_o(0:imax+1,jmax/p_row,kmax/p_col)
      real dft_o(0:imax+1,jmax/p_row,kmax/p_col)
      real dfz_o(0:imax+1,jmax/p_row,kmax/p_col)
      real dfc_o(0:imax+1,jmax/p_row,kmax/p_col)
      real dfr_oo(0:imax+1,jmax/p_row,kmax/p_col)
      real dft_oo(0:imax+1,jmax/p_row,kmax/p_col)
      real dfz_oo(0:imax+1,jmax/p_row,kmax/p_col)
      real dfc_oo(0:imax+1,jmax/p_row,kmax/p_col)
      real fac1,fac2,fac3
      real u_t(0:imx,jmax,kmax/p_col) 
      real v_t(0:imx,jmax,kmax/p_col) 
      real w_t(0:imx,jmax,kmax/p_col) 
      real c_t(0:imx,jmax,kmax/p_col)
      real lambda_t(0:imx,jmax,kmax/p_col)

      real ww1(0:imx,jmax,kmax/p_col) 
      real ww2(0:imx,jmax,kmax/p_col) 
      real ww3(0:imx,jmax,kmax/p_col) 
      real cc(0:imx,jmax,kmax/p_col)
 
      common /ab3/dfr_n,dft_n,dfz_n,dfc_n,dfr_o,dft_o,dfz_o,dfc_o,dfr_oo,dft_oo,dfz_oo,dfc_oo
      integer order
      stime = MPI_WTIME()

      if (ini .eq. 0) then
      call momz(dfr_n,dft_n,dfz_n,dfc_n,ww1,ww2,ww3,cc,unew,vnew,wnew,cnew,
     ^ u_t,v_t,w_t,c_t,visc,lambdanew,lambda_t,imax,jmax,kmax,imx,p_row,p_col,ru,rp,dr,dtheta,dz,nrank,px)
      call momt(dfr_n,dft_n,dfz_n,dfc_n,ww1,ww2,ww3,cc,unew,vnew,wnew,cnew,
     ^ u_t,v_t,w_t,c_t,visc,lambda_t,imax,jmax,kmax,imx,p_row,p_col,ru,rp,dr,dtheta,dz,nrank,px)
      call momr(dfr_n,dft_n,dfz_n,dfc_n,unew,vnew,wnew,cnew,visc,lambdanew,imax,jmax,kmax,imx,
     ^ p_row,p_col,ru,rp,dr,dtheta,dz,mr_s,mr_c,Twb,Twt,nrank,px)
      dfr_o = dfr_n
      dft_o = dft_n
      dfz_o = dfz_n
      dfc_o = dfc_n
      dfr_oo = dfr_n
      dft_oo = dft_n
      dfz_oo = dfz_n
      dfc_oo = dfc_n
      endif 

      if (ini.ne.0) then
      call momz(dfr_n,dft_n,dfz_n,dfc_n,ww1,ww2,ww3,cc,unew,vnew,wnew,cnew,
     ^ u_t,v_t,w_t,c_t,visc,lambdanew,lambda_t,imax,jmax,kmax,imx,p_row,p_col,ru,rp,dr,dtheta,dz,nrank,px)
      call momt(dfr_n,dft_n,dfz_n,dfc_n,ww1,ww2,ww3,cc,unew,vnew,wnew,cnew,
     ^ u_t,v_t,w_t,c_t,visc,lambda_t,imax,jmax,kmax,imx,p_row,p_col,ru,rp,dr,dtheta,dz,nrank,px)
      call momr(dfr_n,dft_n,dfz_n,dfc_n,unew,vnew,wnew,cnew,visc,lambdanew,imax,jmax,kmax,imx,
     ^ p_row,p_col,ru,rp,dr,dtheta,dz,mr_s,mr_c,Twb,Twt,nrank,px)

      endif
      fac1 = 3./2. ! 23./12.
      fac2 =-1./2.  !-4./3.
      fac3 = 0. ! 5./12.
      dudt = unew + dt * (fac1*dfr_n+fac2*dfr_o + fac3*dfr_oo)
      dvdt = vnew + dt * (fac1*dft_n+fac2*dft_o + fac3*dft_oo)
      dwdt = wnew + dt * (fac1*dfz_n+fac2*dfz_o + fac3*dfz_oo   +  dpdx)     
      dcdt = cnew + dt * (fac1*dfc_n+fac2*dfc_o + fac3*dfc_oo + 17.55/(Re*Pr) )


      dfr_oo=dfr_o 
      dft_oo=dft_o 
      dfz_oo=dfz_o
      dfc_oo=dfc_o 

      dfr_o =dfr_n
      dft_o =dft_n
      dfz_o =dfz_n
      dfc_o =dfc_n
c      if (rank.eq.0) then
c         write(6,*)  'CPU time in Adamsb Bashforth3    =  ', MPI_WTIME()-stime 
c      endif
      end



      subroutine Etotexp(time)
      include 'param.txt'
      include 'common.txt'
      real time
      end

      real function ener()
      include 'param.txt'
      include 'common.txt'
      integer istap
      character*5 cha
      real um(imax),vm(imax),wm(imax)
      real ur(imax),vr(imax),wr(imax),uw(imax),ener1
      um=0
      vm=0
      wm=0
      ur=0
      vr=0
      wr=0
      uw=0
      do k=1,kmax/p_col
        do j=1,jmax/p_row
           do i=1,imax
      um(i)=um(i)+unew(i,j,k)/(px*jmax*kmax)
      vm(i)=vm(i)+vnew(i,j,k)/(px*jmax*kmax)
      wm(i)=wm(i)+wnew(i,j,k)/(px*jmax*kmax)
            enddo
         enddo
      enddo
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
            ur(i)=ur(i)+(unew(i,j,k)-um(i))**2
            vr(i)=vr(i)+(vnew(i,j,k)-vm(i))**2
            wr(i)=wr(i)+(wnew(i,j,k)-wm(i))**2
            if (i.lt.imax) uw(i)=uw(i)+unew(i,j,k)*(wnew(i,j,k)+wnew(i+1,j,k)-wm(i)-wm(i+1))*0.5
          enddo
        enddo
      enddo
      ener1 =0
      do i = 1 , imax 
       ener1=ener1+ rp(i)*(ru(i)-ru(i-1))*dtheta*dz*0.5*(ur(i)+vr(i)+wr(i))/(px*jmax*kmax)
      enddo
      ener = ener1
      end

      subroutine cmpbs(bulk,cbulk,str1,str2,fl1,fl2)
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real um(imax),vm(imax),wm(imax),wcm(imax),bulk1,hi_s(imax),hi_c(0:imax),wmmm(imax)
      real ur(imax),vr(imax),wr(imax),uw(imax),cm(imax),ener1,stress1,stress2,flux1,flux2
      real str1,str2,fl1,fl2,bulk,cbulk,cbulk1
      wm=0
      cm=0
      uw=0
      wcm=0
      do k=1,kmax/p_col
      do j=1,jmax/p_row
       do i=1,imax
      wm(i)=wm(i)+wnew(i,j,k)/(jmax*kmax/(p_col*p_row))
      cm(i)=cm(i)+cnew(i,j,k)/(jmax*kmax/(p_col*p_row))
      wcm(i)=wcm(i)+wnew(i,j,k)*cnew(i,j,k)/(jmax*kmax/(p_col*p_row))
            enddo
         enddo
      enddo
      bulk1 = 0
      cbulk1= 0
      do i=1,imax
       bulk1 = bulk1 + (Ru(i)-Ru(i-1))/2*wm(i)
      cbulk1 = cbulk1 + (Ru(i)-Ru(i-1))/2*wcm(i)
      enddo
      call mpi_allreduce(wm,wmmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      wm = wmmm/px
      call mpi_allreduce(cm,wmmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      cm = wmmm/px
      hi_s = cm
      call der1c_s_c6(imax,hi_s,hi_c,dr,Twb,Twt)
      hi_c = hi_c*mr_c
      flux1 = -hi_c(imax)/(Re*Pr)
      flux2 = -hi_c(0)/(Re*Pr)
      hi_s = wm
      call der1w_s_c6(imax,hi_s,hi_c,dr)
      hi_c = hi_c*mr_c
      stress1 = -visc*hi_c(imax)
      stress2 = -visc*hi_c(0)
      call mpi_allreduce(bulk1  ,bulk  ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(cbulk1  ,cbulk,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(stress1,str1  ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(stress2,str2  ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(flux1,fl1     ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(flux2,fl2     ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      bulk = bulk/ px
      cbulk = cbulk/(px*bulk)
      str1 = str1 /px
      str2 = str2 /px
      fl1 = fl1/px
      fl2 = fl2/px
      end

      subroutine chkdt()
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real tmp1,tmp2,dt1,dr1,umax,dtmax
      integer ier
      dtmax =1.5e-4
      umax =0
      do k=1,kmax/p_col
	do j=1,jmax/p_row
	  do i=1,imax 
      dr1 = Ru(i)-Ru(i-1)
      dt1 = dtheta
      
      tmp1 = (abs(unew(i,j,k))/dr1 + abs(vnew(i,j,k))/dt1 + abs(wnew(i,j,k))/dz)
!      tmp2 =10*visc*(1./dr1**2+1./dz**2)
      dt1 = 0.2/(tmp1) 
      umax = max (unew(i,j,k),vnew(i,j,k),wnew(i,j,k),umax)
      enddo
      enddo
      enddo
      call mpi_allreduce(dt1,dt    ,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ier)
      dt  = min ( dt , dtmax)
      end

      subroutine report
      include 'param.txt'
      include 'common.txt'
      write(6,*) '**************************************************************'
      write(6,*) '**************************************************************'
      write(6,*) '***  DNS of turbulent pipe flow                            ***'
      write(6,111) imax,jmax,kmax
 111  format(     '***  resolution (Nr,Nt,Nz) = ' , 3i7    , '            ***')
      write(6,112) px
 112  format(     '***  Number of processors  = ' , i7              , '         ***')
      write(6,*)
      write(6,*) '**************************************************************'
      write(6,*) '***** number of points  = ', imax*jmax*kmax
      end

      function f_colebrook(Re,bulk)
      real cole
      f=0.0001
      do i=1,20
         Reb = Re*bulk
	f = f - cole(f,Reb)/((cole(f+1e-5,Reb)-cole(f-1e-5,Reb))/2e-5)
      enddo
      f_colebrook=f
      end
      function cole(f,Reb)
       cole = ( 1/sqrt(f) +2.*log10(2.51/(Reb*sqrt(f))))
      end
      function f_shokling(Re,bulk)
      real shok 
      f=0.0001
      do i=1,20
         Reb = Re*bulk
	f = f - shok(f,Reb)/((cole(f+1e-5,Reb)-cole(f-1e-5,Reb))/2e-5)
      enddo
      f_shokling=f
      end
      function shok(f,Reb)
       shok= ( 1/sqrt(f) -1.93*log10(Reb*sqrt(f))+0.537)
      end


      subroutine trunc(f,rank)
      use decomp_2d
      include 'param.txt'
      real wt(2*jmax+15),wk(2*kmax+15)
      real f   (0:i1          ,jmax/p_row,kmax/p_col),fz(kmax),ft(jmax),dk(kmax),dj(jmax)
      real t_f (0:imx,jmax    ,kmax/p_col)
      real t_fz(0:imx,jmax/p_col,kmax    )
      integer mask(imax),ifil(0:imax+1),i_dex
      do i=1,imax
         ifil(i)=jmax
      enddo
      idex = i_dex(nrank)
      ifil(0)=ifil(1)
      ifil(imax+1)=ifil(imax)
      do i=1,imax
         if (ifil(i).gt. jmax) ifil(i)=jmax
      enddo
      call vrffti(jmax,wt)
      call vrffti(kmax,wk)
      call transpose_x_to_y(f,t_f)
      do i=0,imx
        do k=1,kmax/p_col
         do j=1,jmax
           ft(j)=t_f(i,j,k)
         enddo
         call vrfftf(1,jmax,ft,dj,1,wt)
         do j=ifil(i+idex),jmax
         ft(j)=0
         enddo
         call vrfftb(1,jmax,ft,dj,1,wt)
          do j=1,jmax
          t_f(i,j,k)=ft(j)
          enddo
        enddo
       enddo
      call transpose_y_to_z(t_f,t_fz)
      do i=0,imx
	do j=1,jmax/p_col
	  do k=1,kmax
	    fz(k)=t_fz(i,j,k)
          enddo
          call vrfftf(1,kmax,fz,dk,1,wk)
          fz(kmax)=0
          call vrfftb(1,kmax,fz,dk,1,wk)
          do k=1,kmax
             t_fz(i,j,k)=fz(k)
          enddo
        enddo
       enddo
       call transpose_z_to_y(t_fz,t_f)
       call transpose_y_to_x(t_f ,  f)
       end

	   

       function i_dex (nrank)
       include 'param.txt'
       i_dex = (nrank/p_col) *(imx+1)
       end  

       function j_dex (nrank)
       include 'param.txt'
       stop
       j_dex = (nrank/p_col) * (jmax/p_row)
       end  

       function k_dex (nrank)
       include 'param.txt'
       stop
       k_dex = nrank/p_row * (kmax/p_col)
       end  
       
      subroutine output2(istap,rank)
      use decomp_2d
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer istap
      character*5 cha
      character*5 cha2
      real um(imax),vm(imax),wm(imax)
      real umm(imax),vmm(imax),wmm(imax)
      real ur(imax),vr(imax),wr(imax),uw(imax),hi_s(imax),hi_c(0:imax)
      real tgrid(0:imax+1,jmax/p_row,kmax/p_col)
      real zgrid(0:imax+1,jmax/p_row,kmax/p_col)
      real t1(0:imx,jmax,kmax/p_col)
      real t2(0:imx,jmax/p_col,kmax)
  

     
      
      integer idex,jdex,kdex,i_dex,j_dex,k_dex
      do k=1,kmax/p_col
	do i=0,imx
	  do j=1,jmax
	    t1(i,j,k)=(J*8*atan(1.)/jmax)
          enddo
        enddo
      enddo
      do k=1,kmax
	do i=0,imx
	  do j=1,jmax/p_col
	    t2(i,j,k)=k*dz
          enddo
        enddo
      enddo
      call transpose_y_to_x(t1,tgrid)
      call transpose_z_to_y(t2,t1)
      call transpose_y_to_x(t1,zgrid)
       
      jdex = j_dex(rank)
      kdex = k_dex(rank)
     
      um=0
      vm=0
      wm=0
      ur=0
      vr=0
      wr=0
      uw=0
      do k=1,kmax/p_col
	 do j=1,jmax/p_row
	    do i=1,imax
	      um(i)=um(i)+unew(i,j,k)/(jmax*kmax*px)
	      vm(i)=vm(i)+vnew(i,j,k)/(jmax*kmax*px)
	      wm(i)=wm(i)+wnew(i,j,k)/(jmax*kmax*px)
	    enddo
	 enddo
      enddo
      call mpi_allreduce(um,umm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(vm,vmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(wm,wmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      um = umm
      vm = vmm
      wm = wmm
      do k=1,kmax/p_col
	do j=1,jmax/p_row
	  do i=1,imax-1
	    ur(i)=ur(i)+(unew(i,j,k)-um(i))**2
	    vr(i)=vr(i)+(vnew(i,j,k)-vm(i))**2
	    wr(i)=wr(i)+(wnew(i,j,k)-wm(i))**2
	    uw(i)=uw(i)+unew(i,j,k)*(wnew(i,j,k)+wnew(i+1,j,k)-wm(i)-wm(i+1))*0.5
	  enddo
	enddo
      enddo
      call mpi_allreduce(ur,umm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(vr,vmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(wr,wmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ur = umm
      vr = vmm
      wr = wmm
      call mpi_allreduce(uw,wmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      uw = wmm 
      hi_s = wm
      call der1w_s_c6(imax,hi_s,hi_c,dr)
      hi_c = hi_c * mr_c
      call cnvstr(istap,cha)
      call cnvstr(rank,cha2)
      open(45,file ='tec_rt3d.'//cha2)
      write(45,*) ' VARIABLES ="X", "Y","Z", "U-vel","V-vel","W-vel","P" '
      write(45,*) ' ZONE I= ',imax, ', J= ', jmax/p_row , 'k = ',kmax/p_col ,'  F=POINT '
      do k=1,kmax/p_col
      do j=1,jmax/p_row
	do i=1,imax
	  write(45,'(6E16.7)') rp(i)*cos(tgrid(i,j,k)),rp(i)*sin(tgrid(i,j,k)),zgrid(i,j,k),
     ^ unew(i,j,k),vnew(i,j,k),wnew(i,j,k),p(i,j,k)
	enddo
      enddo
      enddo
      close(45)
      end
      subroutine wait(n,x)
      x =1000.
      do i=1,n*1000
         x = sin(x)+cos(x)
      enddo
      end
