plot_equations(0.1, 2, 3, 20, 1);
plot_equations(0.5, 10, 3, 20, 2);
plot_equations(0.75, 15, 3, 20, 3);

plot_equations_with_error_norm(0.2, 2, 3, 10);
plot_equations_with_error_norm(0.1, 2, 3, 20);

function [] = plot_equations_with_error_norm(delta_t, t_max, k, num_iter)
N = 0 : 1 : num_iter;
t = 0 : delta_t : t_max;

analytic_sol = exp(-k*t);
close_form_2 = (1-k*delta_t).^N;
close_form_3 = ( 1 + (delta_t / 2) * (-2*k + k^2 * delta_t) ).^N;

error_close_2 = norm(close_form_2 - analytic_sol);
error_close_3 = norm(close_form_3 - analytic_sol);
fprintf('Error_close_2: %.12f \n', error_close_2);
fprintf('Error_close_3: %.12f \n', error_close_3);

% figure(figure_num)
% title('Difference equation analysis')
% plot(t, analytic_sol)
% hold on
% plot(t, close_form_2)
% plot(t, close_form_3)
% ylabel('Function value')
% xlabel('time')
% legend('analytic','closed-form-2', 'closed-form-3')
% hold off

end

function [] = plot_equations(delta_t, t_max, k, num_iter, figure_num)

N = 0 : 1 : num_iter;
t = 0 : delta_t : t_max;

analytic_sol = exp(-k*t);
close_form_2 = (1-k*delta_t).^N;
close_form_3 = ( 1 + (delta_t / 2) * (-2*k + k^2 * delta_t) ).^N;

figure(figure_num)
title('Difference equation analysis')
plot(t, analytic_sol)
hold on
plot(t, close_form_2)
plot(t, close_form_3)
ylabel('Function value')
xlabel('time')
legend('analytic','closed-form-2', 'closed-form-3')
hold off


end