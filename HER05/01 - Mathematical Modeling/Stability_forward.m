function [T_u, T_l, H_u, H_l, Y_u, Y_l, Q_u, Q_l, Mx_u, My_u, Mx_l, My_l, Cq_i_u, Cq_i_l, theta_0_u, theta_0_l] = Stability_forward(u, v, w, rho_input, m_input, g_input, polar)

[T_u, H_u, Y_u, Q_u, Mx_u, My_u, Cq_i_u, theta_0_u] = Stability_forward_upper(u, v, w, rho_input, m_input, g_input, polar);
[T_l, H_l, Y_l, Q_l, Mx_l, My_l, Cq_i_l, theta_0_l] = Stability_forward_lower(u, v, w, rho_input, m_input, g_input, polar);

end