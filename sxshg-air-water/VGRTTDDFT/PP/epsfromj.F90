!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Program epsfromj
  implicit none
  integer :: ihw,ixyz,iter,i,nt0,Nt,Nomega,IOs
  CHARACTER(len=48) :: arg, filen
  CHARACTER*2::ext_field,Longi_Trans
  real(8) :: hw,tt,dt,Eamp(3),domega,avg(1),dAc,pi,tfac,tdamp,ttl,ttr
  complex(8) :: jav_w(1),E_ext_w(1),E_tot_w(1),pav_w(1)
  complex(8) :: zsigma_w(1),zeps(1),zI
  real(8),allocatable ::javt(:,:),tz(:),Pavt(:,:)

  if(iargc() .lt. 2) then
     write(6,*) "Usage: ftz-jw.x file-t Nt [ Nomega domega dAc tdamp ]"
     stop
  endif
  Nomega=2000
  domega=0.001
  tfac=1.0d0
  tdamp=855.d0   ! 855 a.u is ~ 21 fs ~ 0.2 eV broadening
  dAc = 1.0d0
  DO i = 1, iargc()
     CALL getarg(i, arg)
!!$     WRITE (*,*) arg
     if(i .eq. 1) then 
        open(7,file=trim(arg),status="OLD")
        filen=trim(arg)
     endif
     if(i .eq. 2) read( arg, '(i10)' ) Nt
     if(i .eq. 3) read( arg, '(i10)'  ) Nomega
     if(i .eq. 4) read( arg, '(f8.5)' ) domega
     if(i .eq. 5) read( arg, '(f12.8)' ) dAc     
     if(i .eq. 6) read( arg, '(f8.0)' ) tdamp     
  END DO


  ext_field = 'LR'
  Longi_Trans ='Tr'
  zI=DCMPLX(0.0d0,1.0d0)
  pi=3.14159265358979323846264338d0

!  write(*,*) dt, Nt, Eamp, Nomega , domega

  allocate(javt(Nt,1),tz(Nt))
  javt(:,:)=0.0d0

  IOs=0
  Nt=1
  do while (IOs .ge. 0)
     read(7,*,IOSTAT=IOs) tz(Nt), javt(Nt,1) 
     if(IOs .lt. 0) then
        Nt=Nt-1
        exit
     endif
     Nt=Nt+1
  enddo
  
  close(7)

  dt = tz(2)-tz(1)
!!$  write(*,*) Nt, javt(Nt,:)

  write(6,'(A,I8,I8,2f8.5,f8.0)') "Nt, Nomega, domega, dAc, tdamp ",Nt,Nomega, domega, dAc, tdamp

!subtract average
  avg(:)=0
  do iter=2,Nt
     avg(:)=avg(:)+dt*(javt(iter-1,:)+javt(iter,:))/2.0
  enddo
  avg(:)=avg(:)/((Nt-1)*dt)

  do i=1,1
     javt(:,i)=javt(:,i)-avg(i)
  enddo

!!$!damp and smooth
!!$
!!$  do iter=1,Nt
!!$     tt=(iter-1)*dt
!!$     javt(iter,:)= javt(iter,:)*smoothing_t(tt,dt,Nt)
!!$  enddo
!!$  
!!$!subtract average again
!!$  avg(:)=0
!!$  do iter=2,Nt
!!$     avg(:)=avg(:)+dt*(javt(iter-1,:)+javt(iter,:))/2.0
!!$  enddo
!!$  avg(:)=avg(:)/((Nt-1)*dt)
!!$
!!$  do i=1,1
!!$     javt(:,i)=javt(:,i)-avg(i)
!!$  enddo



!  javt = javt/13.60581d0   !convert to a.u?

  filen="ftout_w.dat"
  open(7,file=filen)
  write(7,'(A,f8.5,I8,2f12.6,f8.0)') "# dt, Nomega, domega, dAc, tdamp ",dt,Nomega, domega, dAc, tdamp  
  filen="sigma_w.dat"
  open(8,file=filen)
  write(8,'(A,f8.5,I8,2f12.6,f8.0)') "# dt, Nomega, domega, dAc, tdamp ",dt,Nomega, domega, dAc, tdamp  
  filen="eps_w.dat"
  open(9,file=filen)
  write(9,'(A,f8.5,I8,2f12.6,f8.0)') "# dt, Nomega, domega, dAc, tdamp ",dt,Nomega, domega, dAc, tdamp  
  do ihw=1,Nomega
     hw=(ihw)*domega !+ zI/(tdamp+dt))
     jav_w=0.d0; 
     do iter=2,Nt
        ttl=(iter-2)*dt
        ttr=(iter-1)*dt
        jav_w(:)=jav_w(:)+(javt(iter-1,:)*exp(zI*hw*ttl)*exp(-ttl/tdamp)*smoothing_t(ttl,dt,Nt)+javt(iter,:)*exp(zI*hw*ttr)*exp(-ttr/tdamp)*smoothing_t(ttr,dt,Nt))/2.0 !*smoothing_t(tt,dt,Nt)!*exp(-tt/tdamp)
     enddo
     jav_w(:)=jav_w(:)*dt; 
     if (ext_field == 'LR') then
        if (Longi_Trans == 'Lo')  then 
!!$           zeps(:)=1.d0/(1.d0-E_tot_w(:)/dAc)
        else if (Longi_Trans == 'Tr') then
           zsigma_w(:)=jav_w(:)/dAc
           zeps=1.d0+zI*4.d0*pi*zsigma_w(:)/hw
        end if
     end if


     if (ext_field == 'LR' .and. Longi_Trans == 'Lo') then
!!$           write(7,'(1x,f13.7,6f22.14)') hw&
!!$                &,(real(zeps(ixyz)),ixyz=1,3)&
!!$                &,(imag(zeps(ixyz)),ixyz=1,3)
     else if (ext_field == 'LR' .and. Longi_Trans == 'Tr') then
        write(7,'(1x,f13.7,2f22.14)') hw*27.2114&
             &,(real(jav_w(ixyz)),ixyz=1,1)&
             &,(imag(jav_w(ixyz)),ixyz=1,1)
        write(8,'(1x,f13.7,2f22.14)') hw*27.2114&        
             &,(real(zsigma_w(ixyz)),ixyz=1,1)&
             &,(imag(zsigma_w(ixyz)),ixyz=1,1)
        write(9,'(1x,f13.7,2f22.14)') hw*27.2114&             
             &,(real(zeps(ixyz)),ixyz=1,1)&
             &,(imag(zeps(ixyz)),ixyz=1,1)
     else
!!$           write(7,'(1x,f13.7,18f22.14)') hw&
!!$                &,(real(jav_w(ixyz)),ixyz=1,3)&
!!$                &,(imag(jav_w(ixyz)),ixyz=1,3)&
!!$                &,(real(E_ext_w(ixyz)),ixyz=1,3)&
!!$                &,(imag(E_ext_w(ixyz)),ixyz=1,3)&
!!$                &,(real(E_tot_w(ixyz)),ixyz=1,3)&
!!$                &,(imag(E_tot_w(ixyz)),ixyz=1,3)
     endif
  enddo
  
  close(7)
!!$  filen="pav_t.dat"
!!$  open(8,file=filen)
!!$  write(8,'(A,f8.5,I8)') "# dt, Nt ",dt,Nt
!!$  do iter=1,Nt
!!$     write(8,'(1x,f13.7,f22.14)') (iter-1)*dt,(javt(iter,ixyz),ixyz=1,1)     
!!$  enddo
  close(8)
  close(9)

Contains
  !======
  !======
  Function smoothing_t(tt,dt,Nt)
    implicit none
    real(8),intent(IN) :: tt,dt
    integer,intent(IN) :: Nt
    real(8) :: smoothing_t

!    smoothing_t=1.d0-3.d0*(tt/(Nt*dt))**2+2.d0*(tt/(Nt*dt))**3
    smoothing_t=1.d0-3.d0*(tt/((Nt-1)*dt))**12+2.d0*(tt/((Nt-1)*dt))**18

    return
  End Function smoothing_t
  !======
  !--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
end Program
