% Dados práticos
clear
clc
%load("fotorreceptor_sem_filtro_azul_V1.mat");
%load("fotorreceptor_sem_filtro_azul_V1_CORRIGIDA.mat");
%load("media_R50.mat");
%load("media_R25.mat");
load("Z_LED_corrigida.mat");
%load("fotorreceptor_com_filtro_azul__com_calibracao_V2.mat");


Dtt = Dtt; 
Z = Z_corrigida;
teta = teta_corrigida;

% Dtt = Dtt_50R; %Dtt %Dtt_25R;
% Z = Z_50R; %Z %Z_25R;
% teta = teta_50R; %teta %teta_25R;


% Dtt = Dtt_25R;
% Z = Z_25R;
% teta = teta_25R;

%%
% Seu código original para estimar a função de transferência
% Converter frequência de Hz para rad/s
frequency = 2 * pi * Dtt;

% Converta a fase de graus para radianos
phase = deg2rad(teta); % Se teta_50R já estiver em graus

% Converta impedância (magnitude) e fase para uma resposta complexa
response = Z .* exp(1i * phase);

% Criar objeto de dados de frequência
data = idfrd(response, frequency, 0);

% Especificar a ordem do modelo
np = 2; % número de polos
nz = 2; % número de zeros
delay = 0; % atraso

peso = linspace(10, 1, length(Dtt));

% Crie opções para tfest com pesos ajustados
options = tfestOptions('WeightingFilter', peso');

% Estimar a função de transferência com opções ajustadas
sys = tfest(data, np, nz, delay, options);

% Comparar a resposta do modelo com os dados de entrada
fig = figure;
set(fig, 'Position', [694.6000 145 841.6000 637.6000]);
compare(data, sys);
grid on;
%title('Comparison between Original Data and Electrical Dynamics Transfer Function');
%title('Comparison between Original Data and Optical Dynamics Transfer Function');
title('Comparação entre Dados Originais e Função de Transferência da Dinâmica Elétrica');

% Extrair os coeficientes do numerador e denominador
[num, den] = tfdata(sys, 'v');

% Exibir a função de transferência em formato de função
disp('Função de Transferência Estimada:');
fprintf('Numerador: ');
disp(num);
fprintf('Denominador: ');
disp(den);

% Mostrar a função de transferência em formato de fração de polinômios
sys_tf = tf(num, den);
disp('Função de Transferência no formato:');
sys_tf

% Complemento para obter a função de transferência na forma de multiplicação de fatores
% Obter polos e zeros
[zeros, poles, gain] = zpkdata(sys, 'v');

% Exibir polos, zeros e ganho
disp('Polos: ');
disp(poles);
disp('Zeros: ');
disp(zeros);
disp('Ganho: ');
disp(gain);

% Calcular raízes (frequências naturais) dos polos e zeros
poles_roots = calculateRoots(poles);
zeros_roots = calculateRoots(zeros);

% Exibir raízes (frequências naturais) dos polos e zeros
disp('Raízes dos Polos (rad/s e Hz):');
disp(poles_roots);
disp('Raízes dos Zeros (rad/s e Hz):');
disp(zeros_roots);

% Função para calcular raízes (frequências naturais) a partir de polos e zeros
function roots = calculateRoots(complexValues)
    % Extrair magnitude (frequência natural em rad/s)
    freq_natural_rad_per_sec = abs(complexValues);
    % Converter de rad/s para Hz
    freq_natural_hz = freq_natural_rad_per_sec / (2 * pi);
    roots = table(complexValues, freq_natural_rad_per_sec, freq_natural_hz, ...
                  'VariableNames', {'ComplexValue', 'Frequency_natural_rad_per_sec', 'Frequency_natural_Hz'});
end