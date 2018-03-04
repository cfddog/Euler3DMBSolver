subroutine Splitting_SW(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
use fieldpm
use ctrlpm
use constant
type(BLOCK_TYPE),pointer :: pBlk
integer no_blk
real,pointer,dimension(:,:,:) :: varx,vary,varz
real fnbd(3),fnbd1(3)
real :: epsi=1.e-6
pBlk=>compblock(no_blk)
  if(iinc .eq. 1) then    
     varx=>pblk%xcx; vary=>pblk%xcy; varz=>pblk%xcz 
  elseif(jinc .eq. 1) then
     varx=>pblk%etx; vary=>pblk%ety; varz=>pblk%etz 
  else
     varx=>pblk%ctx; vary=>pblk%cty; varz=>pblk%ctz 
  endif
  
  i=icur
  j=jcur
  k=kcur
  do icnt=nst,ned
     U2=pblk%u(i,j,k)**2+pblk%v(i,j,k)**2+pblk%w(i,j,k)**2
     a=SQRT(gama*pblk%p(i,j,k)/pblk%rho(i,j,k))
     h=gama/(gama-1.)*pblk%p(i,j,k)/pblk%rho(i,j,k)+.5*U2
     Uc=pblk%u(i,j,k)*varx(i,j,k)+pblk%v(i,j,k)*vary(i,j,k)+pblk%w(i,j,k)*varz(i,j,k)
     vmod=SQRT(varx(i,j,k)**2+vary(i,j,k)**2+varz(i,j,k)**2)
     v1=a*vmod
     fnbd(1)=Uc-v1
     fnbd(2)=Uc
     fnbd(3)=Uc+v1
     varx_=varx(i,j,k)/vmod
     vary_=vary(i,j,k)/vmod
     varz_=varz(i,j,k)/vmod
     Uc_=Uc/vmod
     v3=pblk%rho(i,j,k)
     dedrho=pblk%p(i,j,k)/pblk%rho(i,j,k)/(gama-1.)+.5*U2
        do i1=1,3
          fnbd1(i1)=.5*(fnbd(i1)+sqrt(fnbd(i1)**2+epsi))        
        enddo
        v1=(fnbd1(3)+fnbd1(1)-2.*fnbd1(2))*.5/gama
        v2=(fnbd1(3)-fnbd1(1))*.5*a/gama
        fp(1,icnt)=v3*(fnbd1(2)+v1)
        fp(2,icnt)=v3*(fnbd1(2)*pblk%u(i,j,k)+v2*varx_+v1*pblk%u(i,j,k))
        fp(3,icnt)=v3*(fnbd1(2)*pblk%v(i,j,k)+v2*vary_+v1*pblk%v(i,j,k))
        fp(4,icnt)=v3*(fnbd1(2)*pblk%w(i,j,k)+v2*varz_+v1*pblk%w(i,j,k))
        fp(5,icnt)=v3*(fnbd1(2)*dedrho+v2*Uc_+v1*h)
        do i1=1,3
          fnbd1(i1)=fnbd(i1)-fnbd1(i1)        
        enddo
        v1=(fnbd1(3)+fnbd1(1)-2.*fnbd1(2))*.5/gama
        v2=(fnbd1(3)-fnbd1(1))*.5*a/gama
        fn(1,icnt)=v3*(fnbd1(2)+v1)
        fn(2,icnt)=v3*(fnbd1(2)*pblk%u(i,j,k)+v2*varx_+v1*pblk%u(i,j,k))
        fn(3,icnt)=v3*(fnbd1(2)*pblk%v(i,j,k)+v2*vary_+v1*pblk%v(i,j,k))
        fn(4,icnt)=v3*(fnbd1(2)*pblk%w(i,j,k)+v2*varz_+v1*pblk%w(i,j,k))
        fn(5,icnt)=v3*(fnbd1(2)*dedrho+v2*Uc_+v1*h)
        i=i+iinc
        j=j+jinc
        k=k+kinc
   enddo
return
    end subroutine
    
subroutine Splitting_GblMaxEign(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
use fieldpm
use ctrlpm
use constant

integer no_blk
type(BLOCK_TYPE),pointer :: pBlk
real,pointer,dimension(:,:,:) :: varx,vary,varz
real flux(5)
real :: GblMaxEign=0.0

pBlk=>compblock(no_blk)
  if(iinc .eq. 1) then    
     varx=>pblk%xcx; vary=>pblk%xcy; varz=>pblk%xcz 
  elseif(jinc .eq. 1) then
     varx=>pblk%etx; vary=>pblk%ety; varz=>pblk%etz 
  else
     varx=>pblk%ctx; vary=>pblk%cty; varz=>pblk%ctz 
  endif
  
  i=icur;j=jcur;k=kcur
  do icnt=nst,ned
     a=sqrt(gama*pblk%p(i,j,k)/pblk%rho(i,j,k))
     tmpMaxEign=ABS(pblk%u(i,j,k)*varx(i,j,k)+pblk%v(i,j,k)*vary(i,j,k)+&
                    pblk%w(i,j,k)*varz(i,j,k))+a*sqrt(varx(i,j,k)**2+vary(i,j,k)**2+varz(i,j,k)**2)
     GblMaxEign=max(tmpMaxEign,GblMaxEign) 
     i=i+iinc;j=j+jinc;k=k+kinc
  enddo
  
  i=icur;j=jcur;k=kcur
  do icnt=nst,ned
        Uc=pblk%u(i,j,k)*varx(i,j,k)+pblk%v(i,j,k)*vary(i,j,k)+pblk%w(i,j,k)*varz(i,j,k)
        rhoUc=pblk%rho(i,j,k)*Uc
        flux(1)=rhoUc
        flux(2)=(rhoUc*pblk%u(i,j,k)+varx(i,j,k)*pblk%p(i,j,k))
        flux(3)=(rhoUc*pblk%v(i,j,k)+vary(i,j,k)*pblk%p(i,j,k))
        flux(4)=(rhoUc*pblk%w(i,j,k)+varz(i,j,k)*pblk%p(i,j,k))
        flux(5)=(gama/(gama-1)*pblk%p(i,j,k)+.5*pblk%rho(i,j,k)* &
                          (pblk%u(i,j,k)**2+pblk%v(i,j,k)**2+pblk%w(i,j,k)**2))*Uc
        do i1=1,5
          fp(i1,icnt)=.5*(flux(i1)+GblMaxEign*pblk%Q(i1,i,j,k))
          fn(i1,icnt)=flux(i1)-fp(i1,icnt)
        enddo
        i=i+iinc
        j=j+jinc
        k=k+kinc
  enddo
  
return
    end subroutine
 subroutine Splitting_LclMaxEign(nst,ned,icur,jcur,kcur,iinc,jinc,kinc,no_blk)
use fieldpm
use ctrlpm
use constant

integer no_blk
type(BLOCK_TYPE),pointer :: pBlk
real,pointer,dimension(:,:,:) :: varx,vary,varz
real flux(5)
real :: GblMaxEign=0.0

pBlk=>compblock(no_blk)
  if(iinc .eq. 1) then    
     varx=>pblk%xcx; vary=>pblk%xcy; varz=>pblk%xcz 
  elseif(jinc .eq. 1) then
     varx=>pblk%etx; vary=>pblk%ety; varz=>pblk%etz 
  else
     varx=>pblk%ctx; vary=>pblk%cty; varz=>pblk%ctz 
  endif
  
  i=icur;j=jcur;k=kcur
  do icnt=nst,ned
       a=sqrt(gama*pblk%p(i,j,k)/pblk%rho(i,j,k))
       vLclMaxEign=ABS(pblk%u(i,j,k)*varx(i,j,k)+pblk%v(i,j,k)*vary(i,j,k)+&
                      pblk%w(i,j,k)*varz(i,j,k))+a*sqrt(varx(i,j,k)**2+vary(i,j,k)**2+varz(i,j,k)**2)
        Uc=pblk%u(i,j,k)*varx(i,j,k)+pblk%v(i,j,k)*vary(i,j,k)+pblk%w(i,j,k)*varz(i,j,k)
        rhoUc=pblk%rho(i,j,k)*Uc
        flux(1)=rhoUc
        flux(2)=(rhoUc*pblk%u(i,j,k)+varx(i,j,k)*pblk%p(i,j,k))
        flux(3)=(rhoUc*pblk%v(i,j,k)+vary(i,j,k)*pblk%p(i,j,k))
        flux(4)=(rhoUc*pblk%w(i,j,k)+varz(i,j,k)*pblk%p(i,j,k))
        flux(5)=(gama/(gama-1)*pblk%p(i,j,k)+.5*pblk%rho(i,j,k)* &
                          (pblk%u(i,j,k)**2+pblk%v(i,j,k)**2+pblk%w(i,j,k)**2))*Uc
        do i1=1,5
          fp(i1,icnt)=.5*(flux(i1)+vLclMaxEign*pblk%Q(i1,i,j,k))
          fn(i1,icnt)=flux(i1)-fp(i1,icnt)
        enddo
        i=i+iinc
        j=j+jinc
        k=k+kinc
  enddo
  
return
    end subroutine
    
    
  
