function x_l = sobre_amostragem(x,M)
%SOBRE_AMOSTRAGEM funcao que realiza sobre amostragem com fator M
%   x - sinal discreto
%   M - num. inteiro: fator de amostragem
n_amostras = length(x);
x_l = zeros(1, M*n_amostras);
x_l(1:M:end) = x; 
end

