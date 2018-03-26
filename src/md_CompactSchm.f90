module md_CompactSchm
    !Coeffs
    !Kim_Opt4: Jae Wook Kim,AIAA JOURNAL Vol. 41, No. 12, December 2003
     real::ain_1= 0.6511278808920836
     real::ain_2= 0.2487500014377899
     real::ain_3= 0.00614479661269978
     real::bin_1= 0.5775233202590945
     real::cin_1= 0.08953895334666784
    !//
     real::ab0_1=-3.061503488555582
     real::ab0_2= 5.917946021057852
     real::ab0_3= 0.4176795271056629
     real::bb0_1= 5.870156099940824
     real::bb0_2= 3.157271034936285
    !//
     real::ab1_0=-0.5401943305881343
     real::ab1_2= 0.8952361063034303
     real::ab1_3= 0.2553815577627246
     real::ab1_4= 0.007549029394582539
     real::bb1_0= 0.1663921564068434
     real::bb1_2= 0.7162501763222718
     real::bb1_3= 0.08619830787164529
    !//
     real::ab2_0=-0.1327404414078232
     real::ab2_1=-0.6819452549637237
     real::ab2_3= 0.7109139355526556
     real::ab2_4= 0.2459462758541114
     real::ab2_5= 0.003965415751510620
     real::bb2_0= 0.03447751898726934
     real::bb2_1= 0.4406854601950040
     real::bb2_3= 0.6055509079866320
     real::bb2_4= 0.08141498512587530
    contains
   !//Jae Wook Kim,AIAA JOURNAL Vol. 41, No. 12, December 2003
    subroutine CentralCompact_Opt4(nst,ned,varin,varout)
    real varin(nst:ned),varout(nst:ned)
    real arr(nst:ned,5),brr(nst:ned),xrr(nst:ned)
	
    !fulfill the matrix
    !//i=nst
    arr(nst,1)=0.;arr(nst,2)=0.0;arr(nst,3)=1.0;arr(nst,4)=bb0_1;arr(nst,5)=bb0_2
    brr(nst)=ab0_1*(varin(nst+1)-varin(nst))+ab0_2*(varin(nst+2)-varin(nst))+ab0_3*(varin(nst+3)-varin(nst))
    !//i=nst+1
    arr(nst+1,1)=0.;arr(nst+1,2)=bb1_0;arr(nst+1,3)=1.0;arr(nst+1,4)=bb1_2;arr(nst+1,5)=bb1_3
    brr(nst+1)=ab1_0*(varin(nst)-varin(nst+1))+ab1_2*(varin(nst+2)-varin(nst+1))+ab1_3*(varin(nst+3)-varin(nst+1))+ab1_4*(varin(nst+4)-varin(nst+1))
    !//i=nst+2
    arr(nst+2,1)=bb2_0;arr(nst+2,2)=bb2_1;arr(nst+2,3)=1.0;arr(nst+2,4)=bb2_3;arr(nst+2,5)=bb2_4
    brr(nst+2)=ab2_0*(varin(nst)-varin(nst+2))+ab2_1*(varin(nst+1)-varin(nst+2))+ab2_3*(varin(nst+3)-varin(nst+2))+ab2_4*(varin(nst+4)-varin(nst+2))+ab2_5*(varin(nst+5)-varin(nst+2))     
    do icnt=nst+3,ned-3
      arr(icnt,1)=cin_1;arr(icnt,2)=bin_1;arr(icnt,3)=1.0;arr(icnt,4)=bin_1;arr(icnt,5)=cin_1
      brr(icnt)=ain_1*(varin(icnt+1)-varin(icnt-1))+ain_2*(varin(icnt+2)-varin(icnt-2))+ain_3*(varin(icnt+3)-varin(icnt-3))
    enddo
    !//i=ned
    arr(ned,5)=0.;arr(ned,4)=0.0;arr(ned,3)=1.0;arr(ned,2)=bb0_1;arr(ned,1)=bb0_2
    brr(ned)=-ab0_1*(varin(ned-1)-varin(ned))-ab0_2*(varin(ned-2)-varin(ned))-ab0_3*(varin(ned-3)-varin(ned))
    !//i=ned-1
    arr(ned-1,5)=0.;arr(ned-1,4)=bb1_0;arr(ned-1,3)=1.0;arr(ned-1,2)=bb1_2;arr(ned-1,1)=bb1_3
    brr(ned-1)=-ab1_0*(varin(ned)-varin(ned-1))-ab1_2*(varin(ned-2)-varin(ned-1))-ab1_3*(varin(ned-3)-varin(ned-1))-ab1_4*(varin(ned-4)-varin(ned-1))
    !//i=ned-2
    arr(ned-2,5)=bb2_0;arr(ned-2,4)=bb2_1;arr(ned-2,3)=1.0;arr(ned-2,2)=bb2_3;arr(ned-2,1)=bb2_4
    brr(ned-2)=-ab2_0*(varin(ned)-varin(ned-2))-ab2_1*(varin(ned-1)-varin(ned-2))-ab2_3*(varin(ned-3)-varin(ned-2))-ab2_4*(varin(ned-4)-varin(ned-2))-ab2_5*(varin(ned-5)-varin(ned-2))  
    !//call band-eqn solution
    call BandEqnSolve(ned-nst+1,2,2,arr,brr,xrr)
    varout=xrr
	
    return
    end subroutine
	
    end module