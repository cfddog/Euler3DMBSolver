module TYPEDEFINE 
    
type BC_TYPE  
integer              :: no_bc    
integer              :: ibst,ibed,jbst,jbed,kbst,kbed,ibctype
integer              :: ibst_targ,ibed_targ,jbst_targ,jbed_targ,kbst_targ,kbed_targ,no_targblk   !if the bc have a Interface B.C.
end type

type BLOCK_TYPE
integer :: no_blk                                      !block number
integer :: ist,ied,jst,jed,kst,ked                     !contains buffer when multiblock
integer :: icmpst,icmped,jcmpst,jcmped,kcmpst,kcmped   !the actual range
integer :: maxi,maxj,maxk                              !dont contain buffer,the actually max dimensions
!flow field data
real,pointer,dimension (:,:,:)   :: u,v,w,p,rho
real,pointer,dimension (:,:,:,:) :: Q,Qn,dQn,dQ
real,pointer,dimension (:,:,:,:) :: RHS
!grid and its metrics
real,pointer,dimension (:,:,:) :: xcoord,ycoord,zcoord
real,pointer,dimension (:,:,:) :: xcx,xcy,xcz
real,pointer,dimension (:,:,:) :: etx,ety,etz
real,pointer,dimension (:,:,:) :: ctx,cty,ctz
real,pointer,dimension (:,:,:) :: dj
!viscosity 
real,pointer,dimension (:,:,:) :: vnulam,vnutur
!bc conditions
integer                                   :: num_BC
type(BC_TYPE),pointer,dimension(:)        :: BC_MSG
end type

end module
       
module fieldpm
use TYPEDEFINE
 integer numblk,numbc,len_buf,nvar
 type(BLOCK_TYPE),pointer,dimension(:) :: compblock
 !helping vectors and tensors
 integer                               :: MINDIMST,MAXDIMED
 real,pointer,dimension(:,:)           :: fp,fn,hflux,hp,hn
 Real time_begin , time_end1 , time_end2
end module
module freepm
!incoming condition
 real vMainf,Re
 real Rhoinf,Uinf
 real Tinf,Tw
end module

module ctrlpm
!conditions
 real dt,cfl
 integer istart
 integer iSplitSchm
 integer iConSchm
 integer iVisSchm
 integer length_buf
 integer itimeSchm
 integer Nmaxstep,Nsub,iperiod
 integer ihybrid
 integer iproj
 real    eps_solver
end module

module constant
real,parameter :: gama=1.4
real,parameter :: prdt=0.72
real,parameter :: prt=0.9
real,parameter :: kapa=0.41
real,parameter :: eps=1.e-20
end module