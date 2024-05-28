function [x,n] = gera_chirp(t_inicial,t_final,f_final,f_amostragem)
%GERA_CHIRP essa funcao gera um sinal chirp
%   t_inicial    - tempo em segundos para inicio do sinal
%   t_final      - tempo em segundos para fim do sinal
%   f_final      - frq. em Hz para fim do sinal
%   f_amostragem - freq. amostragem

A_0 = 2*pi*(f_final)/(2*(t_final-t_inicial));
T_a = 1/f_amostragem;
n_amostras = (t_final-t_inicial)/T_a;
n = 0:n_amostras-1;
% n = linspace...
x = cos(A_0*(T_a*n).*(T_a*n));
end

