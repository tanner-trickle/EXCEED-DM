* **config**

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

  * :code:`p_vec_list` - List of :math:`\mathbf{p}` vectors for each state.

    * **Dim**: [:math:`N`, 3]
