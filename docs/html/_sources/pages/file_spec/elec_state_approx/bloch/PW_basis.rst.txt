* **config**

  * :code:`G_list_red` - List of :math:`\mathbf{G}` vectors, in reduced coordinates, the Fourier components of the Bloch wave functions, :math:`\widetilde{u}_\mathbf{G}`, are computed for.

    * **Dim**: [:math:`N_\mathbf{G}`, 3]

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

  * **u_FT_c**

    * :code:`n_<n>` - Imaginary part of :math:`\widetilde{u}_\mathbf{G}` for the state :code:`<n>`, :math:`1 \leq \langle n \rangle \leq N`. Spin components can be specified by increasing the second dimension of this dataset to 2, i.e., :math:`N_s = 2`.

      * **Dim**: [:math:`N_\mathbf{G}`, :math:`N_s`]

  * **u_FT_r**

    * :code:`n_<n>` - Real part of :math:`\widetilde{u}_\mathbf{G}` for the state :code:`<n>`, :math:`1 \leq \langle n \rangle \leq N`. Spin components can be specified by increasing the second dimension of this dataset to 2, i.e., :math:`N_s = 2`.

      * **Dim**: [:math:`N_\mathbf{G}`, :math:`N_s`]
