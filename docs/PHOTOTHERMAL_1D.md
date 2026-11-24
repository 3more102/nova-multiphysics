# Planar photothermal reference

The optical Helmholtz cell absorptance A_i is a fraction of incident normal power flux I (W/m²). The cell source is q_i = I A_i / h_i (W/m³), where h_i is cell thickness. For constant conductivity k_i and uniform cell-averaged heating, dF/dx=q_i and F=-k_i dT/dx. The exact cell relations are F_R=F_L+I A_i and T_R=T_L-(h_i/k_i)(F_L+I A_i/2).

The initial face flux follows from the prescribed left and right temperatures. The solver reconstructs nodal temperatures and face heat flux in linear time and reports the absolute signed energy residual F_right-F_left-sum(I A_i). Positive flux points right.

Use the overload accepting Helmholtz1DResult to transfer its cell absorptance without interpolation. The optical and thermal cells must share the same ordered mesh; pass the physical incident flux, not electric-field amplitude. Heating is uniform within each cell, not a reconstruction of the intra-cell optical field. This is a one-way steady planar reference with fixed temperature boundaries, not a 3D thermal FEM solver, transient solver, or self-consistent thermo-optic model.

Tests cover a symmetric heated slab, unequal conductivity, asymmetric heating, the absorption-result adapter, and invalid total absorption. Full repository CI is required before merge.
