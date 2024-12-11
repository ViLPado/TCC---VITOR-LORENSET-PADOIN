function [Z_corrigida,teta_corrigida] = calibracao(Z,teta,Dtt)

% O objetivo desta função é remover o erro do sistema de medição.
% As variáveis de entrada são:
%   Z - vetor do módulo da impedância, gerado pela função "dadoCSVslim"
%   teta - vetor do ângulo da impedância, gerado pela função "dadoCSVslim"
%   Dtt - vetor da frequência da impedância, gerado pela função "dadoCSVslim"
% 
% As variáveis de saída são:
%   Z_corrigida - vetor do módulo da impedância com o erro do sistema
%   removido
%   teta_corrigida - vetor do ângulo da impedância com o erro do sistema
%   removido
% 
% São necessários os arquivos de calibração estarem na pasta do programa em
% que a função "calibrcao" será usada, importo-os com:
% load("media_R50.mat");
% load("media_R25.mat");
% 
% 

%%
%carregar os vetores salvos da média das três aquisições

load("media_R50.mat");
load("media_R25.mat");

%%

% PRA AJUSTRA O TAMANHO DO VETOR DO FATOR D CORREÇÃO PARA TER O MESMO
% TAMNHO DO VETOR DA CARGA TESTADA E CORREPOSNDE AS MESMAS FREQUÊNCIAS
% USADAS PARA COM A CARGA SENDO TESTADA

Dtt = Dtt;

% Pego o comprimento do vetor de referência (que pode ser o Dtt_50R ou Dtt_25R, pois ambos tem 37 pontos)
% e dessa forma posso econtrar o valor mínimo e máximo da frequência da 
% carga que está sendo testada, dessa forma, faço a sdequação do código
% para que o vetor de correção apenas corrija na faixa de frequência
% adequada.
% Por exemplo, o fator de correção corrige entre 1 kHz e 10 MHz. Logo, se
% estou fazendo um teste com uma carga e a faixa de frequência usada é entre 10
% kHz e 1 MHz. Preciso fazer com que o fator de correção tenha o mesmo
% número de pontos da aquisição da carga testada, além de pegar os fator
% para corrigir adequados para cada frequência.


freq_teste_inicial = find(abs(Dtt_50R - min(Dtt)) <= 100);  % Econtro o ponto de frequência min da carga testada que corresponde ao qual localização no vetor do fator de correção, um uma tolerância de +-100Hz
freq_teste_final  = find(abs(Dtt_50R - max(Dtt)) <= 100);   % Econtro o ponto de frequência máx da carga testada que corresponde ao qual localização no vetor do fator de correção, um uma tolerância de +-100Hz

Z_50R = Z_50R(freq_teste_inicial:freq_teste_final);
Z_25R = Z_25R(freq_teste_inicial:freq_teste_final);

teta_50R = teta_50R(freq_teste_inicial:freq_teste_final);
teta_25R = teta_25R(freq_teste_inicial:freq_teste_final);

%%

% PARA CALCULO DO FATOR DE CORREÇÃO NO PONTO MÁXIMO, DE Z USANDO O 50R



Z_fator_mod_max = 50./Z_50R;
Z_fator_deg_max = -teta_50R;

% Convertendo ângulo para radianos
radianos_50R = deg2rad(Z_fator_deg_max);

% Calculando partes real e imaginária
correct_real_50R = Z_fator_mod_max .* cos(radianos_50R);
correct_imaginaria_50R = Z_fator_mod_max .* sin(radianos_50R);

Zfactor_correct_MAX = correct_real_50R + (correct_imaginaria_50R*1i);

%-------------------------------------------------------------------------
% PARA CALCULO DO FATOR DE CORREÇÃO NO PONTO MÁXIMO, DE Z USANDO O 25R

Z_fator_mod_min = 25./Z_25R;
Z_fator_deg_min = -teta_25R;

% Convertendo ângulo para radianos
radianos_25R = deg2rad(Z_fator_deg_min);

% Calculando partes real e imaginária
correct_real_25R = Z_fator_mod_min .* cos(radianos_25R);
correct_imaginaria_25R = Z_fator_mod_min .* sin(radianos_25R);

Zfactor_correct_MIN = correct_real_25R + (correct_imaginaria_25R*1i);



%-----------------------------------------------------------------------
% Convertendo ângulo para radianos

teta = teta;
Z = Z;

radianos = deg2rad(teta);

% Calculando partes real e imaginária
real = Z .* cos(radianos);
imaginaria = Z .* sin(radianos);

Zret = real + (imaginaria*1i);

modulo_valor_medido = Z;

% O valor medido é apenas o múdulo de Z da carga que está sendo testada
proporcao_correcao = (modulo_valor_medido - Z_25R)/(Z_50R - Z_25R);

fator_correcao = (Zfactor_correct_MAX * proporcao_correcao) + (Zfactor_correct_MIN * (1-proporcao_correcao));


Zret = Zret.*fator_correcao %Impedânci corrigida, ou seja, erro do sistema "removido"

Z_corrigida = abs(Zret) %Módulo da impedância corrigida

teta_corrigida = rad2deg(angle(Zret)) %Ângulo da impedância corrigida
