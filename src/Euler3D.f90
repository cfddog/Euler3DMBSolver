!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!                       !
!                 Program for 3D Euler Equations                      ! 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
    program main    
use fieldpm
use freepm
use ctrlpm
  call Readparameters()
  call gridBCIO
  call AllocateVariables()
!//generate the initial field
  call InitialField()
!//compute the metrics
  call ComputeMetrics()
!// primtive parameters to conservation terms
  call PvartoQvar()
!!// call the solvers: euler or NS
  call Solvers()
!!// output results at last
!  call OutputResults(imax,jmax,kmax,u,v,w,p,rho)  
!!// call postprocess
!   call Postprocess(imax,jmax,kmax)
!//end    
end program