    module BCTypeParameters
    !//No BC
    integer,parameter :: itypeBC_NoSpecifyBC =0
    !//Wall BCs
    integer,parameter :: itypeBC_InViscidWall=1
    !//Farfield BCs
    integer,parameter :: itypeBC_FarField    =2
    !//Inlet BCs
    integer,parameter :: itypeBC_SuperInflow =3
    !//Outlet BCs
    integer,parameter :: itypeBC_SuperOutflow=4
    !//Periodical BCs
    integer,parameter :: itypeBC_Period=5
    !//Symmetry BCs
    integer,parameter :: itypeBC_SymmPlane  =6
    integer,parameter :: itypeBC_2DSymmPlane=61
    !//Block-Block interfaces
    integer,parameter :: itypeBC_B2BInterface=-1
    end module
    