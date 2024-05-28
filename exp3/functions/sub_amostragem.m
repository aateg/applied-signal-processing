function x_d = sub_amostragem(x, M)
%SUB_AMOSTRAGEM funcao que realiza sub amostragem por um fator M
%   x - sinal discreto
%   M - num. inteiro: fator de amostragem
n_amostras = length(x);
n_d = 1:M:n_amostras;
x_d = x(n_d);
end

