    program po2d


    implicit none
    integer                 i, j, k
    integer                 Nmesh, Stop_itr
    real                    Re, dt, eps


    write(*,*) 'input meshsize'
    read(*,*) Nmesh

    allocate(psi(Nmesh,Nmesh))
    allocate(omg(Nmesh,Nmesh))
    allocate(tmp(Nmesh,Nmesh))


    do i = 1, Nmesh
      psi(:,:)=0.0
      omg(:,:)=0.0
      tmp(:,:)=0.0
    end do

    write(*,*)'Mesh size=',Nmesh,'*',Nmesh,'. all grid initialized.'
    write(*,*) psi, omg, tmp
    stop
    end
