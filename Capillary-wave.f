program order_parameter
    implicit none
    integer, parameter :: N = 1000  ! Nombre d'atomes (ajuster si nécessaire)
    integer, parameter :: l = 6     ! Degré des harmoniques sphériques
    real(8), parameter :: r_cut = 3.0  ! Rayon de coupure pour les voisins
    real(8), parameter :: pi = 3.141592653589793d0

    ! Déclaration des variables
    real(8) :: positions(N,3), q6(N,13), Sij(N), phi(N)
    integer :: neighbors(N, N), num_neighbors(N)
    integer :: i, j, m
    real(8) :: rij(3), r_norm, theta, phi_angle, dot_product, norm_i, norm_j

    ! Charger les positions depuis un fichier
    open(unit=10, file='positions.xyz', status='old', action='read')
    do i = 1, N
        read(10,*) positions(i,1), positions(i,2), positions(i,3)
    end do
    close(10)

    ! Trouver les voisins
    call find_neighbors(N, positions, neighbors, num_neighbors, r_cut)

    ! Calcul des q6m pour chaque atome
    call compute_q6(N, positions, neighbors, num_neighbors, q6, l)

    ! Calcul de Sij
    call compute_Sij(N, q6, neighbors, num_neighbors, Sij)

    ! Calcul du paramètre d'ordre phi
    call compute_phi(N, Sij, neighbors, num_neighbors, phi)

    ! Écriture des résultats
    open(unit=20, file='order_parameter.dat', status='replace')
    do i = 1, N
        write(20,*) positions(i,3), phi(i)
    end do
    close(20)

    print*, "Calcul terminé. Résultats écrits dans 'order_parameter.dat'"

contains

    ! Trouver les voisins dans un rayon r_cut
    subroutine find_neighbors(N, positions, neighbors, num_neighbors, r_cut)
        implicit none
        integer, intent(in) :: N
        real(8), intent(in) :: positions(N,3)
        integer, intent(out) :: neighbors(N, N), num_neighbors(N)
        real(8), intent(in) :: r_cut
        integer :: i, j
        real(8) :: rij(3), r_norm

        num_neighbors(:) = 0
        neighbors(:,:) = 0

        do i = 1, N
            do j = 1, N
                if (i /= j) then
                    rij = positions(j,:) - positions(i,:)
                    r_norm = sqrt(sum(rij**2))
                    if (r_norm < r_cut) then
                        num_neighbors(i) = num_neighbors(i) + 1
                        neighbors(i, num_neighbors(i)) = j
                    end if
                end if
            end do
        end do
    end subroutine find_neighbors

    ! Calcul des q6m
    subroutine compute_q6(N, positions, neighbors, num_neighbors, q6, l)
        implicit none
        integer, intent(in) :: N, l
        real(8), intent(in) :: positions(N,3)
        integer, intent(in) :: neighbors(N, N), num_neighbors(N)
        real(8), intent(out) :: q6(N, 13)
%Lire les positions des atomes.
%Trouver les voisins de chaque atome avec un critère de distance.
%Calculer les harmoniques sphériques q6(i)q6​(i).
%Calculer le facteur de corrélation SijSij​.
%Calculer le paramètre d’ordre ϕiϕi​.
%Écrire les résultats dans un fichier.


        integer :: i, j, m
        real(8) :: rij(3), r_norm, theta, phi_angle

        q6(:,:) = 0.0d0

        do i = 1, N
            if (num_neighbors(i) == 0) cycle

            do j = 1, num_neighbors(i)
                rij = positions(neighbors(i,j),:) - positions(i,:)
                r_norm = sqrt(sum(rij**2))
                theta = acos(rij(3) / r_norm)
                phi_angle = atan2(rij(2), rij(1))

                do m = -6, 6
                    q6(i, m+7) = q6(i, m+7) + sph_harm(m, l, phi_angle, theta)
                end do
            end do

            q6(i,:) = q6(i,:) / num_neighbors(i)
        end do
    end subroutine compute_q6

    ! Fonction pour les harmoniques sphériques (simplifiée pour Fortran)
    function sph_harm(m, l, phi, theta) result(value)
        implicit none
        integer, intent(in) :: m, l
        real(8), intent(in) :: phi, theta
        real(8) :: value
        value = cos(m * phi) * sin(l * theta)  ! Approximation basique
    end function sph_harm

    ! Calcul de Sij
    subroutine compute_Sij(N, q6, neighbors, num_neighbors, Sij)
        implicit none
        integer, intent(in) :: N
        real(8), intent(in) :: q6(N, 13)
        integer, intent(in) :: neighbors(N, N), num_neighbors(N)
        real(8), intent(out) :: Sij(N)
        integer :: i, j
        real(8) :: dot_product, norm_i, norm_j

        Sij(:) = 0.0d0

        do i = 1, N
            if (num_neighbors(i) == 0) cycle

            do j = 1, num_neighbors(i)
                dot_product = sum(q6(i,:) * q6(neighbors(i,j),:))
                norm_i = sqrt(sum(q6(i,:)**2))
                norm_j = sqrt(sum(q6(neighbors(i,j),:)**2))

                if (norm_i > 0.0d0 .and. norm_j > 0.0d0) then
                    Sij(i) = Sij(i) + dot_product / (norm_i * norm_j)
                end if
            end do

            Sij(i) = Sij(i) / num_neighbors(i)
        end do
    end subroutine compute_Sij

    ! Calcul du paramètre d'ordre φ_i
    subroutine compute_phi(N, Sij, neighbors, num_neighbors, phi)
        implicit none
        integer, intent(in) :: N
        real(8), intent(in) :: Sij(N)
        integer, intent(in) :: neighbors(N, N), num_neighbors(N)
        real(8), intent(out) :: phi(N)
        integer :: i, j

        phi(:) = 0.0d0

        do i = 1, N
            if (num_neighbors(i) == 0) cycle

            do j = 1, num_neighbors(i)
                phi(i) = phi(i) + Sij(neighbors(i,j))
            end do

            phi(i) = phi(i) / num_neighbors(i)
        end do
    end subroutine compute_phi

end program order_parameter
