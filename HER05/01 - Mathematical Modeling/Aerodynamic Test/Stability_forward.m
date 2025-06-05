function [T_u, T_l, H_u, H_l, Y_u, Y_l, Q_u, Q_l, Mx_u, My_u, Mx_l, My_l, Cq_i_u, Cq_i_l, beta1c_UR, beta1s_UR, beta1c_LR, beta1s_LR] = Stability_forward(u, v, w, a_input, rho_input, mu_input, m_input, g_input, polar, theta_0_u, theta_0_l, theta_1c, theta_1s)

[T_u, H_u, Y_u, Q_u, Mx_u, My_u, Cq_i_u, beta1c_UR, beta1s_UR] = Stability_forward_upper(u, v, w, a_input, rho_input, mu_input, m_input, g_input, polar, theta_0_u, theta_0_l, theta_1c, theta_1s);
[T_l, H_l, Y_l, Q_l, Mx_l, My_l, Cq_i_l, beta1c_LR, beta1s_LR] = Stability_forward_lower(u, v, w, a_input, rho_input, mu_input, m_input, g_input, polar, theta_0_u, theta_0_l, theta_1c, theta_1s);

end