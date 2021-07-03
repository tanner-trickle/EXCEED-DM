module transition
    !! Handles labelling of the transitions that each processor should compute 
    !! for

    implicit none

    integer :: n_tran
        !! Total number of transitions
        !!
        !! = n_init*n_fin

    integer :: n_tran_per_proc
        !! Number of transitiosn each processor has to calculate

    integer, allocatable :: tran_to_init_fin_id(:, :)
        !! Dim : [n_tran, 2]
        !!
        !! Each transition (i, f) is given a unique index, and this 
        !! is the map back. For each transition id return the initial 
        !! or final state index  

    integer, allocatable :: job_table(:, :)
        !! Dim : [n_proc, n_tran_per_proc]
        !! give each processor a list of transitions to compute for

contains

    subroutine set_job_table(n_proc, n_init, n_fin, verbose)
        implicit none

        integer :: n_proc

        integer :: n_init, n_fin

        logical, optional :: verbose

        integer :: i, j, tran_id, id, f

        if ( verbose ) then

            print*, 'Configuring jobs for processors...'
            print*

        end if

        n_tran = n_init*n_fin

        allocate(tran_to_init_fin_id(n_tran, 2))

        id = 0
        do i = 1, n_init
            do f = 1, n_fin

                id = id + 1

                tran_to_init_fin_id(id, 1) = i
                tran_to_init_fin_id(id, 2) = f

            end do
        end do 

        if ( mod(n_tran, n_proc) .eq. 0 ) then 
            n_tran_per_proc = n_tran/n_proc
        else
            n_tran_per_proc = n_tran/n_proc + 1
        end if 

        allocate(job_table(n_proc, n_tran_per_proc))

        tran_id = 0
        
        do j = 1, n_tran_per_proc
            do i = 1, n_proc
                
                tran_id = tran_id + 1

                if (tran_id .gt. n_tran) then 
                    job_table(i, j) = 0
                else
                    job_table(i, j) = tran_id
                end if 

            end do
        end do  
        
        if ( verbose ) then             

            print*, '----------------------------------------'
            print*
            
            if ( mod(n_tran, n_proc) .eq. 0 ) then 
        
                print*, '    Equal processor load.'
                print*
            
            else if ( n_tran_per_proc .eq. 1 ) then 

                print*, '    Number of processors is greater than the number of i -> f',&
                        '    transitions. Consider lowering the number of processors.'
                print*
                print*, '    Number of transitions = ', n_tran
                print*

            else

                print*, '    Unequal processor load. Some processors will be given null jobs.'
                print*
                print*, '    Number of transitions = ', n_tran
                print*

            end if 

            print*, '    Number of calculations per processor = ', n_tran_per_proc
            print*   
            print*, '----------------------------------------'
            print*
            
        end if

    end subroutine

end module
