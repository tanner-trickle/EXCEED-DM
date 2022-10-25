* **config**

  * :code:`n_r_vec_grid` - Number of :math:`\mathbf{r}` vectors to keep in each direction. For example, :code:`n_r_vec_grid = [3, 3, 3]` keeps all :math:`\mathbf{r} \in \{ r_1, r, 2, r_3\}, -1 \leq r_i \leq 1`.

    * **Dim**: [3]

  * :code:`n_x_grid` - Number of :math:`\mathbf{x}` points in the unit cell to compute the wave functions for (in each direction).  :math:`\begin{align} \mathbf{x}_{ijk} = \{ \frac{i - 1}{N_1}, \frac{j - 1}{N_2}, \frac{k - 1}{N_3} \} \end{align}` in reduced coordinates.

    * **Dim**: [3]

* **state_info**

  * :code:`Zeff_list` - List of :math:`Z_\text{eff}` parameters for each state. Only used when transitioning to a :code:`bloch/single_PW` state. 

    * **Dim**: [:math:`N`]

  * :code:`energy_list` - Energy, :math:`\omega`, of each state.

    * **Units**: :math:`\text{eV}`

    * **Dim**: [:math:`N`]

  * :code:`i_list` - Band index, :math:`i`, of each state.

    * **Dim**: [:math:`N`]

  * :code:`jac_list` - List of jacobians, :math:`j_I`, of discretized state sum. Let :math:`N'` be the exact number of initial states. Therefore, in theory, to sum over the states one applies :math:`(1/N') \sum_{I'}`. However in practice :math:`N'` is usually very, very large, and the sum must be discretized, :math:`(1/N') \sum_{I'} \rightarrow \sum_I j_I` and this array contains that jacobian for each state.

    * **Dim**: [:math:`N`]

  * :code:`k_id_list` - The :math:`\mathbf{k}` index of each state, i.e., :math:`i` if the :math:`\mathbf{k}` vectors are indexed by :math:`i`.

    * **Dim**: [:math:`N`]

  * :code:`k_vec_red_list` - :math:`\mathbf{k}` vectors for each state in reduced coordinates.

    * **Dim**: [:math:`N`, 3]

  * :code:`coeff_list` - Collection of STO coefficients for each state. Each state is represented by 4 length :math:`N_j` vectors, :math:`\{ \mathbf{n}, \mathbf{Z}, \mathbf{N}, \mathbf{C} \}`. The :math:`\mathbf{n}` coefficients for the nth state are stored in [n, :, 1]. The :math:`\mathbf{Z}` coefficients for the nth state are stored in [n, :, 2]. The :math:`\mathbf{N}` coefficients for the nth state are stored in [n, :, 3]. The :math:`\mathbf{C}` coefficients for the nth state are stored in [n, :, 4].

    * **Dim**: [:math:`N`, :math:`N_j`, 4]

  * :code:`eq_pos_red_list` - Equilibrium position (in reduced coordinates) of the ion each state is associated with.

    * **Dim**: [:math:`N`, 3]

  * :code:`nj_list` - Number of expansion coefficients for each state. 

    * **Dim**: [:math:`N`]

  * :code:`nlm_list` - :math:`n, \ell, m` quantum numbers of each state, stored in that order.

    * **Dim**: [:math:`N`, 3]
