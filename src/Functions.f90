subroutine PvarToQvar()
use constant
use fieldpm
type(BLOCK_TYPE),pointer :: pBlk
do iblk=1,numblk
    pBlk=>compblock(iblk)
  do i=pBlk%ist,pBlk%ied
    do j=pBlk%jst,pBlk%jed
      do k=pBlk%kst,pBlk%ked
        pBlk%Q(1,i,j,k)=pBlk%rho(i,j,k)
        pBlk%Q(2,i,j,k)=pBlk%rho(i,j,k)*pBlk%u(i,j,k)
        pBlk%Q(3,i,j,k)=pBlk%rho(i,j,k)*pBlk%v(i,j,k)
        pBlk%Q(4,i,j,k)=pBlk%rho(i,j,k)*pBlk%w(i,j,k)
        pBlk%Q(5,i,j,k)=pBlk%p(i,j,k)/(gama-1)+0.5*pBlk%rho(i,j,k)*(pBlk%u(i,j,k)**2+pBlk%v(i,j,k)**2+pBlk%w(i,j,k)**2)
        do i1=1,5
            pBlk%Qn(i1,i,j,k)=pBlk%Q(i1,i,j,k)
        enddo

      enddo
     enddo
  enddo
  
enddo

return
end subroutine  
!
subroutine QvarToPvar()
use constant
use fieldpm
type(BLOCK_TYPE),pointer :: pBlk
do iblk=1,numblk
    pBlk=>compblock(iblk)
  do i=pBlk%ist,pBlk%ied
    do j=pBlk%jst,pBlk%jed
      do k=pBlk%kst,pBlk%ked
        pBlk%rho(i,j,k)=pBlk%Q(1,i,j,k)
        pBlk%u(i,j,k)=pBlk%Q(2,i,j,k)/pBlk%Q(1,i,j,k)
        pBlk%v(i,j,k)=pBlk%Q(3,i,j,k)/pBlk%Q(1,i,j,k)
        pBlk%w(i,j,k)=pBlk%Q(4,i,j,k)/pBlk%Q(1,i,j,k) 
        pBlk%p(i,j,k)=(gama-1)*( pBlk%Q(5,i,j,k)-0.5*( pBlk%Q(2,i,j,k)**2+pBlk%Q(3,i,j,k)**2+pBlk%Q(4,i,j,k)**2 )/pBlk%Q(1,i,j,k) )      
        enddo
     enddo
   enddo
 enddo
return
end subroutine 
 subroutine Readparameters()
 use fieldpm
 use ctrlpm
 use freepm
 namelist /FreeParamters/ vMainf,Re,rhoinf,uinf,Tinf,Tw
 namelist /CtrlParamters/ isolver,nvar,dt,cfl,iSplitSchm,iConSchm,IVisSchm,iTimeSchm,&
                          Nmaxstep,Nsub,iperiod,ihybrid,iproj,eps_solver,istart,len_buf
 open(10,file="INPUT.in")
 read(10,*) !freestream conditions
 read(10,*) vMainf,Re,Tinf,Tw
 read(10,*) !Equation sets
 read(10,*) isolver,nvar
 read(10,*) !controls
 read(10,*) iSplitSchm,iConSchm,IVisSchm,iTimeSchm
 read(10,*) dt,Nmaxstep,iperiod
 read(10,*) istart
 read(10,*) eps_solver
 close(10)
 
 return  
 end subroutine
 subroutine AllocateVariables()
 use fieldpm
 type(BLOCK_TYPE),pointer :: pBlk
 !reading the grid.get the dimensions
 do iblk=1,numblk
     pblk=>compblock(iblk)
     allocate( pBlk%u(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked),&
               pBlk%v(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked),&
               pBlk%w(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked),&
               pBlk%p(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked),&
               pBlk%rho(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked) )
			   
     allocate( pBlk%Q(nvar,pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked),&
               pBlk%Qn(nvar,pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked) )
     allocate( pBlk%dQn(nvar,pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked),&
               pBlk%dQ(nvar,pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked) )
			   
     allocate( pBlk%RHSi(nvar,pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked) )
	 allocate( pBlk%RHSj(nvar,pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked) )
	 allocate( pBlk%RHSk(nvar,pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked) )
	 allocate( pBlk%RHS(nvar,pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked) )
     
     allocate( pBlk%xcx(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked),&
               pBlk%xcy(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked),&
               pBlk%xcz(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked) )
     allocate( pBlk%etx(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked),&
               pBlk%ety(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked),&
               pBlk%etz(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked) ) 
     allocate( pBlk%ctx(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked),&
               pBlk%cty(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked),&
               pBlk%ctz(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked) ) 
     allocate( pBlk%dj(pBlk%ist:pBlk%ied,pBlk%jst:pBlk%jed,pBlk%kst:pBlk%ked) )
 enddo
 !for global variables
    allocate(fp(nvar,MINDIMST:MAXDIMED))
    allocate(fn(nvar,MINDIMST:MAXDIMED))
    allocate(hflux(nvar,MINDIMST-1:MAXDIMED))
    allocate(hp(nvar,MINDIMST:MAXDIMED))
    allocate(hn(nvar,MINDIMST:MAXDIMED))
    allocate(flux(nvar,MINDIMST:MAXDIMED))
 return
    end subroutine 

subroutine InitialField()
use fieldpm
use Freepm
use constant
type(BLOCK_TYPE),pointer :: pBlk
!
!  if(istart .eq. 1) then                   !reading the restart file
     
  if(istart .eq. 2) then             !reading the existed field file as the start field 
    open(10,file="output_RES.dat",form="formatted")
       read(10,*) nb
       do iblk=1,numblk
           pBlk=>compblock(iblk)
           read(10,*) id,jd,kd
       enddo
       do iblk=1,numblk
           pBlk=>compblock(iblk)
           read(10,*) d1,d2,d3,d4
           read(10,*) (((pBlk%rho(i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
           (((pBlk%u(i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
           (((pBlk%v(i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
           (((pBlk%w(i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
           (((pBlk%p(i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped)
       enddo
       close(10)   
  else  !istart =0                              !init from the zero initialfield
    do iblk=1,numblk
       pBlk=>compblock(iblk)
       do i=pBlk%ist,pBlk%ied
          do j=pBlk%jst,pBlk%jed
              do k=pBlk%kst,pBlk%ked 
                 if(pblk%xcoord(i,j,k) .le. 0) then
                   pBlk%rho(i,j,k)=0.125
                   pBlk%u(i,j,k)=0.0
                   pBlk%v(i,j,k)=0.0
                   pBlk%w(i,j,k)=0.0
                   pBlk%p(i,j,k)=0.1
                else
                   pBlk%rho(i,j,k)=1.0
                   pBlk%u(i,j,k)=0.0
                   pBlk%v(i,j,k)=0.0
                   pBlk%w(i,j,k)=0.0
                   pBlk%p(i,j,k)=1.0                 
                endif
              enddo
          enddo
       enddo  
    enddo    
 endif    
return          
    end subroutine   
    
 subroutine OutputResults()
 use fieldpm
 
 type(BLOCK_TYPE),pointer :: pBlk
!the plot3D flow field
 open(10,file="output.dat",form="formatted")
    write(10,*) numblk
    do iblk=1,numblk
        pBlk=>compblock(iblk)
        write(10,*) pblk%maxi,pblk%maxj,pblk%maxk
    enddo
    do iblk=1,numblk
        pBlk=>compblock(iblk)
        write(10,*) 0.0,0.0,0.0,0.0
        write(10,*) (((pBlk%rho(i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
        (((pBlk%u(i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
        (((pBlk%v(i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
        (((pBlk%w(i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped),&
        (((pBlk%p(i,j,k),i=pBlk%icmpst,pblk%icmped),j=pblk%jcmpst,pblk%jcmped),k=pblk%kcmpst,pblk%kcmped)
    enddo
    close(10)

    !open(10,file="output_RHS.dat",form="formatted")
    !write(10,*) numblk
    !do iblk=1,numblk
    !    pBlk=>compblock(iblk)
    !    write(10,*) pblk%ied-pblk%ist+1,pblk%jed-pblk%jst+1,pblk%ked-pblk%kst+1
    !enddo
    !do iblk=1,numblk
    !    pBlk=>compblock(iblk)
    !    write(10,*) 0.0,0.0,0.0,0.0
    !    write(10,*) (((pBlk%RHS(1,i,j,k),i=pBlk%ist,pblk%ied),j=pblk%jst,pblk%jed),k=pblk%kst,pblk%ked),&
    !    (((pBlk%RHS(2,i,j,k),i=pBlk%ist,pblk%ied),j=pblk%jst,pblk%jed),k=pblk%kst,pblk%ked),&
    !    (((pBlk%RHS(3,i,j,k),i=pBlk%ist,pblk%ied),j=pblk%jst,pblk%jed),k=pblk%kst,pblk%ked),&
    !    (((pBlk%RHS(4,i,j,k),i=pBlk%ist,pblk%ied),j=pblk%jst,pblk%jed),k=pblk%kst,pblk%ked),&
    !    (((pBlk%RHS(5,i,j,k),i=pBlk%ist,pblk%ied),j=pblk%jst,pblk%jed),k=pblk%kst,pblk%ked)
    !enddo
    !close(10)
    
    end subroutine
    !//
    subroutine RLmaxtrix(no_blk,icur,jcur,kcur,iinc,jinc,kinc,RM,LM)
    use fieldpm
    use constant
    type(BLOCK_TYPE),pointer :: pBlk
    real,pointer,dimension(:,:,:) :: varx,vary,varz
    integer no_blk
    real LM(5,5),RM(5,5)
    pblk=>compblock(no_blk)
    
    if(iinc .eq. 1) then    
        varx=>pblk%xcx; vary=>pblk%xcy; varz=>pblk%xcz 
    elseif(jinc .eq. 1) then
        varx=>pblk%etx; vary=>pblk%ety; varz=>pblk%etz 
    else
        varx=>pblk%ctx; vary=>pblk%cty; varz=>pblk%ctz 
    endif
    i=icur;j=jcur;k=kcur
    !//Roe averaging
    iL=icur
    jL=jcur
    kL=kcur
    !//
    iR=icur+iinc
    jR=jcur+jinc
    kR=kcur+kinc
    !//
    rhoL=pblk%rho(iL,jL,kL)
    rhoR=pblk%rho(iR,jR,kR)
    rhoc=((sqrt(rhoL)+sqrt(rhoR))/2.)**2
    !//
    uL=pblk%u(iL,jL,kL)
    uR=pblk%u(iR,jR,kR)
    uc=(sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR))
    !//
    vL=pblk%v(iL,jL,kL)
    vR=pblk%v(iR,jR,kR)
    vc=(sqrt(rhoL)*vL+sqrt(rhoR)*vR)/(sqrt(rhoL)+sqrt(rhoR))
    !//
    wL=pblk%w(iL,jL,kL)
    wR=pblk%w(iR,jR,kR)
    wc=(sqrt(rhoL)*wL+sqrt(rhoR)*wR)/(sqrt(rhoL)+sqrt(rhoR))
    !//
    pL=pblk%p(iL,jL,kL)
    pR=pblk%p(iR,jR,kR)
    HL=gama*pL/(gama-1)/rhoL+0.5*(uL**2+vL**2+wL**2)
    HR=gama*pR/(gama-1)/rhoR+0.5*(uR**2+vR**2+wR**2)
    Hc=(sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR))
    !//
    pc=(gama-1)/gama*(rhoc*Hc-0.5*rhoc*(uc**2+vc**2+wc**2))
    cc=sqrt((gama-1)*(Hc-0.5*(uc**2+vc**2+wc**2)))
    !//
     vkx=(sqrt(rhoL)*varx(iL,jL,kL)+sqrt(rhoR)*varx(iR,jR,kR))/(sqrt(rhoL)+sqrt(rhoR))
     vky=(sqrt(rhoL)*vary(iL,jL,kL)+sqrt(rhoR)*vary(iR,jR,kR))/(sqrt(rhoL)+sqrt(rhoR))
     vkz=(sqrt(rhoL)*varz(iL,jL,kL)+sqrt(rhoR)*varz(iR,jR,kR))/(sqrt(rhoL)+sqrt(rhoR))
    !//L and R characteristic matrix
     vkmod=sqrt(vkx**2+vky**2+vkz**2)
     vkx_=vkx/vkmod
     vky_=vky/vkmod
     vkz_=vkz/vkmod
     sita=vkx*uc+vky*vc+vkz*wc
     sita_=sita/vkmod
     Xi=0.5*(uc**2+vc**2+wc**2)/cc
     phi2=0.5*(gama-1)*(uc**2+vc**2+wc**2)
     !
     LM(1,1)= vkx_/cc*(cc*cc-phi2)+vkz_*vc-vky_*wc
     LM(1,2)= vkx_/cc*(gama-1)*uc
     LM(1,3)= vkx_/cc*(gama-1)*vc-vkz_
     LM(1,4)= vkx_/cc*(gama-1)*wc+vky_
     LM(1,5)=-vkx_/cc*(gama-1)
     !
     LM(2,1)= vky_/cc*(cc*cc-phi2)+vkx_*wc-vkz_*uc
     LM(2,2)= vky_/cc*(gama-1)*uc+vkz_
     LM(2,3)= vky_/cc*(gama-1)*vc
     LM(2,4)= vky_/cc*(gama-1)*wc-vkx_
     LM(2,5)=-vky_/cc*(gama-1)
     !
     LM(3,1)= vkz_/cc*(cc*cc-phi2)+vky_*uc-vkx_*vc
     LM(3,2)= vkz_/cc*(gama-1)*uc-vky_
     LM(3,3)= vkz_/cc*(gama-1)*vc+vkx_
     LM(3,4)= vkz_/cc*(gama-1)*wc
     LM(3,5)=-vkz_/cc*(gama-1)
     !
     LM(4,1)= 1/sqrt(2.)/cc*(phi2-cc*sita_)
     LM(4,2)=-1/sqrt(2.)/cc*((gama-1)*uc-vkx_*cc)
     LM(4,3)=-1/sqrt(2.)/cc*((gama-1)*vc-vky_*cc)
     LM(4,4)=-1/sqrt(2.)/cc*((gama-1)*wc-vkz_*cc)
     LM(4,5)= 1/sqrt(2.)/cc*(gama-1)
     !
     LM(5,1)= 1/sqrt(2.)/cc*(phi2+cc*sita_)
     LM(5,2)=-1/sqrt(2.)/cc*((gama-1)*uc+vkx_*cc)
     LM(5,3)=-1/sqrt(2.)/cc*((gama-1)*vc+vky_*cc)
     LM(5,4)=-1/sqrt(2.)/cc*((gama-1)*wc+vkz_*cc)
     LM(5,5)= 1/sqrt(2.)/cc*(gama-1)
     !//
     RM(1,1)= vkx_/cc
     RM(1,2)= vky_/cc
     RM(1,3)= vkz_/cc
     RM(1,4)= 1./sqrt(2.)/cc
     RM(1,5)= 1./sqrt(2.)/cc
     !
     RM(2,1)= vkx_/cc*uc
     RM(2,2)= vky_/cc*uc+vkz_
     RM(2,3)= vkz_/cc*uc-vky_
     RM(2,4)= 1./sqrt(2.)/cc*(uc+vkx_*cc)
     RM(2,5)= 1./sqrt(2.)/cc*(uc-vkx_*cc)
     !
     RM(3,1)= vkx_/cc*vc-vkz_
     RM(3,2)= vky_/cc*vc
     RM(3,3)= vkz_/cc*vc+vkx_
     RM(3,4)= 1./sqrt(2.)/cc*(vc+vky_*cc)
     RM(3,5)= 1./sqrt(2.)/cc*(vc-vky_*cc)
     !
     RM(4,1)= vkx_/cc*wc+vky_
     RM(4,2)= vky_/cc*wc-vkx_
     RM(4,3)= vkz_/cc*wc
     RM(4,4)= 1./sqrt(2.)/cc*(wc+vkz_*cc)
     RM(4,5)= 1./sqrt(2.)/cc*(wc-vkz_*cc)
     !
     RM(5,1)= Xi*vkx_+wc*vky_-vc*vkz_
     RM(5,2)= Xi*vky_+uc*vkz_-wc*vkx_
     RM(5,3)= Xi*vkz_+vc*vkx_-uc*vky_
     RM(5,4)= 1./sqrt(2.)/cc*((phi2+cc*cc)/(gama-1)+sita_*cc)
     RM(5,5)= 1./sqrt(2.)/cc*((phi2+cc*cc)/(gama-1)-sita_*cc)
     !!//test
     do i1=1,5
         do i2=1,5
             tmp=0.0
             do ii=1,5
                 tmp=tmp+LM(i1,ii)*RM(ii,i2)
             enddo
             if(i1 .eq. i2) then
                 if(tmp .gt. 1.+5.e-5 .or. tmp .le. 1.- 5.e-5) then
                     print*,i1,i2,tmp
                     pause
                 endif
             endif
         enddo
     enddo
     
    return
    end subroutine
