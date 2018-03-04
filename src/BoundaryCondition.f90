!this will be a complex problem,but for now,we just use it simply
subroutine EnforceBC()
use constant
use fieldpm
use freepm
use BCTypeParameters
type(BLOCK_TYPE),pointer :: pblk
type(BC_TYPE),pointer :: pbc

do iblk=1,numblk
    pblk=>compblock(iblk)
   do ibc=1,pblk%num_BC
       pbc=>pblk%BC_MSG(ibc)
       select case (pbc%ibctype)
       case(itypeBC_NoSpecifyBC)
           
       case(itypeBC_InViscidWall)
            call BCProcess_InvisidWall(iblk,pbc%ibst,pbc%ibed,pbc%jbst,pbc%jbed,pbc%kbst,pbc%kbed)
       case(itypeBC_FarField)
           call BCProcess_FarField(iblk,pbc%ibst,pbc%ibed,pbc%jbst,pbc%jbed,pbc%kbst,pbc%kbed)
       case(itypeBC_SuperInflow)
           call BCProcess_SuperInflow(iblk,pbc%ibst,pbc%ibed,pbc%jbst,pbc%jbed,pbc%kbst,pbc%kbed)
       case(itypeBC_SuperOutflow)
           call BCProcess_SuperOutflow(iblk,pbc%ibst,pbc%ibed,pbc%jbst,pbc%jbed,pbc%kbst,pbc%kbed)
       case(itypeBC_Period)
       case(itypeBC_2DSymmPlane)
           call BCProcess_2DSymmPlane(iblk,pbc%ibst,pbc%ibed,pbc%jbst,pbc%jbed,pbc%kbst,pbc%kbed)
       case(itypeBC_B2BInterface)
           call BCProcess_B2BInterface(pbc%ibst,pbc%ibed,pbc%jbst,pbc%jbed,pbc%kbst,pbc%kbed,iblk,&
                                       pbc%ibst_targ,pbc%ibed_targ,pbc%jbst_targ,pbc%jbed_targ,pbc%kbst_targ,&
                                       pbc%kbed_targ,pbc%no_targblk)
       case default
          !//noBC used as default
       end select
   enddo 
enddo

return
    end subroutine
    
    subroutine BCProcess_SuperInflow(noblk,ibst,ibed,jbst,jbed,kbst,kbed)
    use fieldpm
    use freepm
    use constant
    type(BLOCK_TYPE),pointer :: pblk
     pblk=>compblock(noblk)
     if(ibst .eq. ibed) then
         iinc=1;jinc=0;kinc=0
     elseif(jbst .eq. jbed) then
         iinc=0;jinc=1;kinc=0
     else
         iinc=0;jinc=0;kinc=1
     endif
     
     do i=ibst,ibed
         do j=jbst,jbed
             do k=kbst,kbed
                 do ivar=1,nvar
                    pBlk%Q(ivar,i,j,k)=pBlk%Qn(ivar,i,j,k)
                 enddo
                 pBlk%rho(i,j,k)=pBlk%Q(1,i,j,k)
                 pBlk%u(i,j,k)  =pBlk%Q(2,i,j,k)/pBlk%Q(1,i,j,k)
                 pBlk%v(i,j,k)  =pBlk%Q(3,i,j,k)/pBlk%Q(1,i,j,k)
                 pBlk%w(i,j,k)  =pBlk%Q(4,i,j,k)/pBlk%Q(1,i,j,k) 
                 pBlk%p(i,j,k)  =(gama-1)*( pBlk%Q(5,i,j,k)-0.5*( pBlk%Q(2,i,j,k)**2+pBlk%Q(3,i,j,k)**2+pBlk%Q(4,i,j,k)**2 )/pBlk%Q(1,i,j,k) )
             enddo
         enddo
     enddo
     
     return
    end subroutine
    !//
    subroutine BCProcess_SuperOutflow(noblk,ibst,ibed,jbst,jbed,kbst,kbed)
    use fieldpm
    use freepm
    use constant
    type(BLOCK_TYPE),pointer :: pblk
     pblk=>compblock(noblk)
 !judge the direction    
     if(ibst .eq. ibed) then
         iinc=1;jinc=0;kinc=0
     elseif(jbst .eq. jbed) then
         iinc=0;jinc=1;kinc=0
     else
         iinc=0;jinc=0;kinc=1
     endif
 !obtained the incre
     if(abs(ibst) .eq. 1) then
         isgn=1
     else
         isgn=-1
     endif
     if(abs(jbst) .eq. 1) then
         jsgn=1
     else
         jsgn=-1
     endif
     if(abs(kbst) .eq. 1) then
         ksgn=1
     else
         ksgn=-1
     endif
     !..
     do i=abs(ibst),abs(ibed),isgn
         do j=abs(jbst),abs(jbed),jsgn
             do k=abs(kbst),abs(kbed),ksgn
                 do ivar=1,nvar
                    pBlk%Q(ivar,i,j,k)=2*pBlk%Q(ivar,i+iinc*isgn,j+jinc*jsgn,k+kinc*ksgn)-pBlk%Q(ivar,i+2*iinc*isgn,j+2*jinc*jsgn,k+2*kinc*ksgn)
                 enddo
                 pBlk%rho(i,j,k)=pBlk%Q(1,i,j,k)
                 pBlk%u(i,j,k)  =pBlk%Q(2,i,j,k)/pBlk%Q(1,i,j,k)
                 pBlk%v(i,j,k)  =pBlk%Q(3,i,j,k)/pBlk%Q(1,i,j,k)
                 pBlk%w(i,j,k)  =pBlk%Q(4,i,j,k)/pBlk%Q(1,i,j,k) 
                 pBlk%p(i,j,k)  =(gama-1)*( pBlk%Q(5,i,j,k)-0.5*( pBlk%Q(2,i,j,k)**2+pBlk%Q(3,i,j,k)**2+pBlk%Q(4,i,j,k)**2 )/pBlk%Q(1,i,j,k) )
             enddo
         enddo
     enddo
     
     return
    end subroutine
    !//
    subroutine BCProcess_InvisidWall(noblk,ibst,ibed,jbst,jbed,kbst,kbed)
    use fieldpm
    use freepm
    use constant
    type(BLOCK_TYPE),pointer :: pblk
    real,pointer,dimension(:,:,:) :: varx,vary,varz
     pblk=>compblock(noblk)
 !judge the direction    
     if(ibst .eq. ibed) then
         iinc=1;jinc=0;kinc=0
     elseif(jbst .eq. jbed) then
         iinc=0;jinc=1;kinc=0
     else
         iinc=0;jinc=0;kinc=1
     endif
 !obtained the incre
     if(abs(ibst) .eq. 1) then
         isgn=1
     else
         isgn=-1
     endif
     if(abs(jbst) .eq. 1) then
         jsgn=1
     else
         jsgn=-1
     endif
     if(abs(kbst) .eq. 1) then
         ksgn=1
     else
         ksgn=-1
     endif
     !/
     if(iinc .eq. 1) then    
     varx=>pblk%xcx; vary=>pblk%xcy; varz=>pblk%xcz 
     elseif(jinc .eq. 1) then
     varx=>pblk%etx; vary=>pblk%ety; varz=>pblk%etz 
     else
     varx=>pblk%ctx; vary=>pblk%cty; varz=>pblk%ctz 
     endif
     do i=abs(ibst),abs(ibed),isgn
         do j=abs(jbst),abs(jbed),jsgn
             do k=abs(kbst),abs(kbed),ksgn
                 sig=sqrt(varx(i,j,k)**2+vary(i,j,k)**2+varz(i,j,k)**2)
                 vnx_=varx(i,j,k)/sig
                 vny_=vary(i,j,k)/sig
                 vnz_=varz(i,j,k)/sig
                 Un=pblk%u(i,j,k)*vnx_+pblk%v(i,j,k)*vny_+pblk%w(i,j,k)*vnz_
                 pblk%u(i,j,k)=pblk%u(i,j,k)-Un*vnx_
                 pblk%v(i,j,k)=pblk%v(i,j,k)-Un*vny_
                 pblk%w(i,j,k)=pblk%w(i,j,k)-Un*vnz_
                 !//
                 pblk%p(i,j,k)=(4.*pblk%p(i+iinc*isgn,j+jinc*jsgn,k+kinc*ksgn)-pblk%p(i+2*iinc*isgn,j+2*jinc*jsgn,k+2*kinc*ksgn))/3.
                 pblk%rho(i,j,k)=(4.*pblk%rho(i+iinc*isgn,j+jinc*jsgn,k+kinc*ksgn)-pblk%rho(i+2*iinc*isgn,j+2*jinc*jsgn,k+2*kinc*ksgn))/3.
                 !//
                 pBlk%Q(1,i,j,k)=pBlk%rho(i,j,k)
                 pBlk%Q(2,i,j,k)=pBlk%rho(i,j,k)*pBlk%u(i,j,k)
                 pBlk%Q(3,i,j,k)=pBlk%rho(i,j,k)*pBlk%v(i,j,k)
                 pBlk%Q(4,i,j,k)=pBlk%rho(i,j,k)*pBlk%w(i,j,k)
                 pBlk%Q(5,i,j,k)=pBlk%p(i,j,k)/(gama-1)+0.5*pBlk%rho(i,j,k)*(pBlk%u(i,j,k)**2+pBlk%v(i,j,k)**2+pBlk%w(i,j,k)**2)
             enddo
         enddo
     enddo
     
     return
    end subroutine
    !//
    subroutine BCProcess_FarField(noblk,ibst,ibed,jbst,jbed,kbst,kbed)
    use fieldpm
    use freepm
    use constant
    type(BLOCK_TYPE),pointer :: pblk
     pblk=>compblock(noblk)
     if(ibst .eq. ibed) then
         iinc=1;jinc=0;kinc=0
     elseif(jbst .eq. jbed) then
         iinc=0;jinc=1;kinc=0
     else
         iinc=0;jinc=0;kinc=1
     endif
     
     do i=ibst,ibed
         do j=jbst,jbed
             do k=kbst,kbed
                 do ivar=1,nvar
                    pBlk%Q(ivar,i,j,k)=pBlk%Qn(ivar,i,j,k)
                 enddo
                 pBlk%rho(i,j,k)=pBlk%Q(1,i,j,k)
                 pBlk%u(i,j,k)  =pBlk%Q(2,i,j,k)/pBlk%Q(1,i,j,k)
                 pBlk%v(i,j,k)  =pBlk%Q(3,i,j,k)/pBlk%Q(1,i,j,k)
                 pBlk%w(i,j,k)  =pBlk%Q(4,i,j,k)/pBlk%Q(1,i,j,k) 
                 pBlk%p(i,j,k)  =(gama-1)*( pBlk%Q(5,i,j,k)-0.5*( pBlk%Q(2,i,j,k)**2+pBlk%Q(3,i,j,k)**2+pBlk%Q(4,i,j,k)**2 )/pBlk%Q(1,i,j,k) )
             enddo
         enddo
     enddo
     
     return
    end subroutine
    
    subroutine BCProcess_2DSymmPlane(noblk,ibst,ibed,jbst,jbed,kbst,kbed)
    use fieldpm
    use freepm
    use constant
    !supposed that there is only THREE LAYER in this direction,take the second plane as the value
    type(BLOCK_TYPE),pointer :: pblk
     pblk=>compblock(noblk)
     if(ibst .eq. ibed) then
         iinc=1;jinc=0;kinc=0
     elseif(jbst .eq. jbed) then
         iinc=0;jinc=1;kinc=0
     else
         iinc=0;jinc=0;kinc=1
     endif
     imid=(pBlk%ied-pBlk%ist)/2+1
     jmid=(pBlk%jed-pBlk%jst)/2+1
     kmid=(pBlk%ked-pBlk%kst)/2+1
     do i=ibst,ibed
         do j=jbst,jbed
             do k=kbst,kbed
                     ii=i*(1-iinc)+iinc*imid
                     jj=j*(1-jinc)+jinc*jmid
                     kk=k*(1-kinc)+kinc*kmid
                 do ivar=1,nvar
                     pBlk%Q(ivar,i,j,k)=pBlk%Q(ivar,ii,jj,kk)
                 enddo
                 pBlk%rho(i,j,k)=pBlk%Q(1,i,j,k)
                 pBlk%u(i,j,k)  =pBlk%Q(2,i,j,k)/pBlk%Q(1,i,j,k)
                 pBlk%v(i,j,k)  =pBlk%Q(3,i,j,k)/pBlk%Q(1,i,j,k)
                 pBlk%w(i,j,k)  =pBlk%Q(4,i,j,k)/pBlk%Q(1,i,j,k) 
                 pBlk%p(i,j,k)  =(gama-1)*( pBlk%Q(5,i,j,k)-0.5*&
                                 ( pBlk%Q(2,i,j,k)**2+pBlk%Q(3,i,j,k)**2+pBlk%Q(4,i,j,k)**2 )/pBlk%Q(1,i,j,k) )
             enddo
         enddo
     enddo
     
     return
     end subroutine
                 
    
    