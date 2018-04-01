module md_CompactSchm
    !Coeffs
	!Kim_OSOT4: Jae Wook Kim,Journal of Computational Acoustics, Vol. 5, No. 2 (1997) 177-191
	 real::Ain_1= 1.568098211519709
     real::Ain_2= 0.2716571074522698
     real::Ain_3=-0.02257678073547548
     real::BETA_IN= 0.4085892691182515
	 !//AT BCP=1
	 real::AB1_1=-2.673444389108146
     real::AB1_2= 1.468066764967325
     real::AB1_3= 1.382688702485047
	 REAL::AB1_4=-0.1773110783442254
	 real::BETA_B1_2=2.701510934904742
	 !//AT BCP=2
	 real::AB2_1=-0.5088675754573845
     real::AB2_2=-0.7029878533366753
     real::AB2_3= 1.040385365448375
	 REAL::AB2_4= 0.1867472036506759
	 REAL::AB2_5=-0.01527714030499072
	 real::BETA_B2_1=0.1532048781838751
	 real::BETA_B2_3=0.7237110491082636
	 !//AT BCP=3
	 real::AB3_1=-0.013127263621621
     real::AB3_2=-0.6038029221734134
     real::AB3_3=-0.4395154246847092
	 REAL::AB3_4= 0.96090920472974
	 REAL::AB3_5= 0.1010303485585628
	 REAL::AB3_6=-0.005493942808558833
	 real::BETA_B3_2=0.2234544771621557
	 real::BETA_B3_4=0.5530910456756884
    contains
   !//Jae Wook Kim,AIAA JOURNAL Vol. 41, No. 12, December 2003
    subroutine CentralCompact_Opt4(nst,ned,varin,varout)
    real varin(nst:ned),varout(nst:ned)
    real arr(nst:ned,3),brr(nst:ned),xrr(nst:ned)
	
    !fulfill the matrix
    !//i=nst
    arr(nst,1)=0.;         arr(nst,2)=1.0;  arr(nst,3)=BETA_B1_2
	brr(nst,1)=AB1_1*varin(nst)+AB1_2*varin(nst+1)+AB1_3*varin(nst+2)+AB1_4*varin(nst+3)
    !//i=nst+1
    arr(nst+1,1)=BETA_B2_1;arr(nst+1,2)=1.0;arr(nst+1,3)=BETA_B2_3
	brr(nst+1,1)=AB2_1*varin(nst)+AB2_2*varin(nst+1)+AB2_3*varin(nst+2)+AB2_4*varin(nst+3)+AB2_5*varin(nst+4)
    !//i=nst+2
    arr(nst+2,1)=BETA_B3_2;arr(nst+2,2)=1.0;arr(nst+2,3)=BETA_B3_4
	brr(nst+2,1)=AB3_1*varin(nst)+AB3_2*varin(nst+1)+AB3_3*varin(nst+2)+AB3_4*varin(nst+3)+AB3_5*varin(nst+4)+AB3_6*varin(nst+5)
    do icnt=nst+3,ned-3
      arr(icnt,1)=BETA_IN;arr(icnt,2)=1.0;arr(icnt,3)=BETA_IN
	  brr(icnt,1)=Ain_1*(varin(icnt+1)-varin(icnt-1))/3.+Ain_2*(varin(icnt+2)-varin(icnt-2))/4.+Ain_3*(varin(icnt+3)-varin(icnt-3))/6.
    enddo
    !//i=ned
    arr(ned,3)=0.0 ;arr(ned,2)=1.0;arr(ned,1)=-BETA_B1_2
	brr(ned,1)=-AB1_1*varin(ned)-AB1_2*varin(ned-1)-AB1_3*varin(ned-2)-AB1_4*varin(ned-3)
    !//i=ned-1
    arr(ned-1,3)=BETA_B2_1;arr(ned-1,2)=1.0;arr(ned-1,1)=BETA_B2_3
	brr(ned-1,1)=-AB2_1*varin(ned)-AB2_2*varin(ned-1)-AB2_3*varin(ned-2)-AB2_4*varin(ned-3)-AB2_5*varin(ned-4)
    !//i=ned-2
    arr(ned-2,3)=BETA_B3_2;arr(ned-2,2)=1.0;arr(ned-2,1)=BETA_B3_4
	brr(ned-2,1)=-AB3_1*varin(ned)-AB3_2*varin(ned-1)-AB3_3*varin(ned-2)-AB3_4*varin(ned-3)-AB3_5*varin(ned-4)-AB3_6*varin(ned-5)
    !//call band-eqn solution
    call BandEqnSolve(ned-nst+1,1,1,arr,brr,xrr)
    varout=xrr
	
    return
    end subroutine
    
    subroutine Explicit_C4(nst,ned,varin,varout)
    real varin(nst:ned),varout(nst:ned)
    do icnt=nst,ned
        if(icnt .eq. nst) then
            varout(icnt)=-11./6.*varin(icnt)+3.*varin(icnt+1)-3./2.*varin(icnt+2)+1./3.*varin(icnt+3)
        elseif(icnt .eq. nst+1) then
            varout(icnt)=-1./3.*varin(icnt-1)-1./2.*varin(icnt)+varin(icnt+1)-1./6*varin(icnt+2)
        elseif(icnt .eq. ned) then
            varout(icnt)=11./6.*varin(icnt)-3.*varin(icnt-1)+3./2.*varin(icnt-2)-1./3.*varin(icnt-3)
        elseif(icnt .eq. ned-1) then
           varout(icnt)= 1./3.*varin(icnt+1)+1./2.*varin(icnt)-varin(icnt-1)+1./6*varin(icnt-2)
        else
            varout(icnt)= 1./12.*varin(icnt-2.)-2./3.*varin(icnt-1)+2./3.*varin(icnt+1)-1./12*varin(icnt+2)
        endif    
    enddo
    
    return
    endsubroutine
    end module