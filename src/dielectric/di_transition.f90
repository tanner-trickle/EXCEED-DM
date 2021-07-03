module di_transition
    !! Handles labelling of the transitions that each processor should compute 
    !! for to compute the dielectric.
    !!
    !! TODO: Make a better interface to transition such that the job table
    !! is not specific to the module. This will allow the same transition
    !! subroutines to be available more generally.

    implicit none

    integer :: di_n_tran
        !! Total number of transitions
        !!
        !! = n_init*n_fin

    integer :: di_n_tran_per_proc
        !! Number of transitiosn each processor has to calculate

    integer, allocatable :: di_tran_to_init_fin_id(:, :)
        !! Dim : [n_tran, 2]
        !!
        !! Each transition (i, f) is given a unique index, and this 
        !! is the map back. For each transition id return the initial 
        !! or final state index  

    integer, allocatable :: di_job_table(:, :)
        !! Dim : [n_proc, n_tran_per_proc]
        !! give each processor a list of transitions to compute for

contains

    subroutine di_set_job_table(n_proc, n_init, n_fin, verbose)
        implicit none

        integer :: n_proc

        integer :: n_init, n_fin

        logical, optional :: verbose

        integer :: i, j, tran_id, id, f

        if ( verbose ) then

            print*, 'Configuring jobs for processors...'
            print*

        end if

        di_n_tran = n_init*n_fin

        allocate(di_tran_to_init_fin_id(di_n_tran, 2))

        id = 0
        do i = 1, n_init
            do f = 1, n_fin

                id = id + 1

                di_tran_to_init_fin_id(id, 1) = i
                di_tran_to_init_fin_id(id, 2) = f

            end do
        end do 

        if ( mod(di_n_tran, n_proc) .eq. 0 ) then 
            di_n_tran_per_proc = di_n_tran/n_proc
        else
            di_n_tran_per_proc = di_n_tran/n_proc + 1
        end if 

        allocate(di_job_table(n_proc, di_n_tran_per_proc))

        tran_id = 0
        
        do j = 1, di_n_tran_per_proc
            do i = 1, n_proc
                
                tran_id = tran_id + 1

                if (tran_id .gt. di_n_tran) then 
                    di_job_table(i, j) = 0
                else
                    di_job_table(i, j) = tran_id
                end if 

            end do
        end do  
        
        if ( verbose ) then             

            print*, '----------------------------------------'
            print*
            
            if ( mod(di_n_tran, n_proc) .eq. 0 ) then 
        
                print*, '    Equal processor load.'
                print*
            
            else if ( di_n_tran_per_proc .eq. 1 ) then 

                print*, '    Number of processors is greater than the number of i -> f',&
                        '    transitions. Consider lowering the number of processors.'
                print*
                print*, '    Number of transitions = ', di_n_tran
                print*

            else

                print*, '    Unequal processor load. Some processors will be given null jobs.'
                print*
                print*, '    Number of transitions = ', di_n_tran
                print*

            end if 

            print*, '    Number of calculations per processor = ', di_n_tran_per_proc
            print*   
            print*, '----------------------------------------'
            print*
            
        end if

    end subroutine

end module
