clear;
clc;

%% 4.1
[sinal, fs] = audioread('antarctica.wav'); % wavread for old versions MATLAB
trecho = sinal(200:439, 1);

[ak, sig2]=lpc(trecho.*hamming(240), 10);

%% 4.2
G = 1;
freqz(G, ak, 512);
hold on;
periodogram(trecho,[],512); % 4.3

%% 4.4
[ak, sig2]=lpc(trecho.*hamming(240), 80); % varie  o valor aqui 20, 40 ou 80
G = 1;
freqz(1, G*ak, 512);
hold on;
periodogram(trecho,[],512);

%% 4.5
plot_fig = false;
pitch = yaapt(sinal,fs,1,[],plot_fig,1);
pitch = [0 pitch];

%% 4.6
N = length(sinal);
sinal_gerado = [];
N_amostras = 80;
N_trecho = 280;
i = 0;
for start=1:N_amostras:(N-N_trecho+N_amostras)
    % 6.a
    i = i + 1;
    if start + N_trecho - 1 > length(sinal) % start+280 > 11201
        trecho = sinal(start:end, 1);
        [ak, sig2] = lpc(trecho.*hamming(length(trecho)), 10);
    else
        trecho = sinal(start:(start+N_trecho-1), 1);
        [ak, sig2] = lpc(trecho.*hamming(N_trecho), 10);
    end

    % 6.b
    if pitch(i) > 0
        G = sqrt(pitch(i)*sig2); % 6.d
        %1. gerar uma sequencia de pulsos
        T = 1/pitch(i); % s
        s = zeros(N_amostras, 1); %0:10ms
        s(1) = 1;
        fator = 10e-3/N_amostras; 
        soma = 0;
        for j=1:N_amostras
            if soma >= T
                s(j) = 1;
                soma = 0;
            end
            soma = soma + fator;
        end
    else
        G = sqrt(sig2); % 6.d
        %2. gere um trecho de ruido com 80 amostras usando randn
        s = randn(N_amostras, 1);
    end
    
    trechosint = filter(G, ak, s);
    sinal_gerado = [sinal_gerado; trechosint];
end

figure;
plot(sinal_gerado);

%% 4.7 Com parâmetros do quantizados
N = length(sinal);
sinal_gerado_quant = [];
N_amostras = 80;
N_trecho = 280;
i = 0;
Ba = 7;
Bg = 5;
for start=1:N_amostras:(N-N_trecho+N_amostras)
    % 6.a
    i = i + 1;
    if start + N_trecho - 1 > length(sinal) % start+280 > 11201
        trecho = sinal(start:end, 1);
        [ak, sig2] = lpc(trecho.*hamming(length(trecho)), 10);
    else
        trecho = sinal(start:(start+N_trecho-1), 1);
        [ak, sig2] = lpc(trecho.*hamming(N_trecho), 10);
    end
    ak = quantize3(ak,Ba);

    % 6.b
    Tp = quantize3(pitch(i), Bg);
    if Tp > 0
        G = quantize3(sqrt(Tp*sig2), Bg);
        
        %1. gerar uma sequencia de pulsos
        fp = 1/Tp; % s - freq pitch
        s = zeros(N_amostras, 1); %0:10ms
        s(1) = 1;
        fator = 10e-3/N_amostras; 
        soma = 0;
        for j=1:N_amostras
            if soma >= fp
                s(j) = 1;
                soma = 0;
            end
            soma = soma + fator;
        end
    else
        G = quantize3(sqrt(sig2), Bg);
        %2. gere um trecho de ruido com 80 amostras usando randn
        s = randn(N_amostras, 1);
    end
    
    trechosint = filter(G, ak, s);
    sinal_gerado_quant = [sinal_gerado_quant; trechosint];
end

hold on; grid on;
plot(sinal_gerado_quant);

%% Parte 6 - Realimentacao
lambda = 1;
[sinal, fs] = audioread('antarctica.wav');

% definicao de variaveis
fa = 8000; % Hz
Ltrecho = 240; % num de amostras dos quadros (30ms)
Npasso = 80; % num de amostras para trechos (10ms)
Q = randn(Npasso, 512); % 512 funcoes aleatorias
K = 2;  % num de funcoes base
p = 10; % ordem do LPC

zs = zeros(p);
sinal_gerado = [];

% processando os quadros
for start=1:Npasso:length(sinal)
    trecho = sinal(start:(start+Ltrecho-1), 1);
    [a_trecho, ~] = lpc(trecho.*hamming(Ltrecho), p);

    % subquadro de comprimento 80 
    % usando as amostras centrais do quadro atual
    subquadro = trecho((Npasso+1):(2*Npasso+1), 1);
    
    % c) filtre todas as sequencias de Q por ak
    Q_filtrado = filter(1, a_trecho, Q);

    % d) melhore a transicao entre os quadros
    [y0, zs] = filter(1, aq, zeros(Npasso,1), zs);

    % e) calcule o sinal e_0
    e_0 = subquadro - y0;

    % f) 
    [ganhos, indices] = find_Nbest_components(e0, Q_filtrado, K); % sinal, codebook_vectors, N

    % g) defina o sinal de excitacao
    d = Q_filtrado(indices,:) * ganhos;

    % h)
    [sinal_saida, zs] = filter(1, ai, d, zs);
    sinal_gerado = [sinal_gerado,; sinal_saida];
end

plot(sinal_gerado);
%% 6-1
fa = 8000; % Hz
L = 240; % num de amostras dos quadros (30ms)
N = 80; % num de amostras para trechos (10ms)
%% 6-2
Q = randn(N, 512); % 512 funcoes aleatorias
K = 2;  % num de funcoes base

%% 6-3
zs = zeros(10); % QUAL O TAMANHO? % condicao inicial do filtro de trato vocal
% acho q eh 10 por cause do p; certo
%% 6-4

% a)
p = 10; 
[aq, sig2] = lpc(quadros(1,:).*hamming(N), p);

% b)
% defina sub-quadro

% c)
[yq, ~] = filter(1, aq, Q); % FIXME

% d) melhorando a transicao entre os quadros
[y0, zs] = filter(1, aq, zeros(N,1), zs);

% e)
e0 = subquadro - y0;

% f)
[ganhos, indices] = find_Nbest_components(e0, Q, K); % sinal, codebook_vectors, N

% g) defina o sinal de excitacao
d = Q(indices,:) * ganhos;

% h)
[sinal_saida, zs] = filter(1, ai, d, zs);

% i) 
beta_quant = quantize3(ganhos, 5); % quantizado com 5bits
ak_quant = quantize3(ak, 7); % quantize o preditor com 7bits

d = Q(indices,:) * beta_quant;

[sinal_saida_quant, zs] = filter(1, ak_quant, d, zs);