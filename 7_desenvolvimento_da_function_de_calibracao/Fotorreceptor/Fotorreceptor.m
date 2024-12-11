%Teste fotorreceptor

tic;
clc
clear


% Listar os arquivos da pasta de interesse

Arquivo = "D:\dados_usuario\Documents\GEDRE\TCC_Git\Código Matlab\7_desenvolvimento_da_function_de_calibracao\Fotorreceptor\Fotorreceptor\Com_flitro_azul\V2";

APENAS_ARQUIVOS_CSV = "\*.csv";

y = dir(Arquivo + APENAS_ARQUIVOS_CSV);

Pot = zeros(1,length(y));
S = zeros(1,length(y));
FP = zeros(1,length(y));
teta = zeros(1,length(y));
Z = zeros(1,length(y));
Zret = zeros(1,length(y));
freq_fft = zeros(1,length(y));
freq_t = zeros(1,length(y));
Dtt = zeros(1,length(y));

% Exibir o caminho de cada arquivo .txt na pasta
for i = 1:length(y)
    caminho = [y(i).folder '\' y(i).name];
    disp(y(i).name);
    [Pot(i),S(i),FP(i),teta(i),Z(i),Zret(i),freq_fft(i),freq_t(i), Dtt(i)] = dadoCSVslim(caminho,1);
end

%Apenas para não ficar rodando o programa toda vez, irei já salvar os
%vetores das aquisições
% OBS: NÃO SÃO OS VALORES CORRIGIDOS, SÃO ELES "BRUTOS"

%save("fotorreceptor_sem_filtro_azul_V1", "Pot","S","FP","teta","Z","Zret","freq_fft","freq_t", "Dtt");
%load("fotorreceptor_sem_filtro_azul_V1");

Amostra = 1:length(y);


tab = table(Amostra',Pot', S', FP', teta', Z', Zret', freq_fft', freq_t', Dtt', 'VariableNames', {'Amostras','Pot', 'S', 'FP', 'teta', 'Z', 'Zret', 'freq pela fft', 'freq pelo delta tempo', 'freq pelo todo tempo'});

tab_mean = table({'Média'}, mean(Pot), mean(S), mean(FP), mean(teta), mean(Z), mean(Zret), mean(freq_fft), mean(freq_t), mean(Dtt), 'VariableNames', {'Amostras','Pot', 'S', 'FP', 'teta', 'Z', 'Zret', 'freq pela fft', 'freq pelo delta tempo', 'freq pelo todo tempo'});

disp(tab);
disp(tab_mean);


% Plotar os dados
% Criar uma figura
figure;

% Primeiro subplot
h(1) = subplot(2,1,1);
semilogx(Dtt, Z, 'o-');
xlabel('Frequência (Hz)');
ylabel('Módulo Z (\Omega)');
%ylim([49 max(Z)]);
title('Módulo Z vs Frequência');
grid on;

% Definir os rótulos do eixo x para mostrar os valores desejados
xticks([10000, 100000, 1000000, 10000000]); % Define os pontos do eixo x
xticklabels({'10kHz', '100kHz', '1MHz', '10MHz'}); % Define os rótulos correspondentes

% Ajustar a aparência dos rótulos do eixo x
xtickangle(45); % Rotaciona os rótulos para melhor legibilidade, se necessário

% Segundo subplot
h(2) = subplot(2,1,2);
semilogx(Dtt, teta, 'o-');
xlabel('Frequência (Hz)');
ylabel('Ângulo teta (°)');
%ylim([-10 max(teta)]);
title('Ângulo teta vs Frequência');
grid on;

% Definir os rótulos do eixo x para mostrar os valores desejados
xticks([10000, 100000, 1000000, 10000000]); % Define os pontos do eixo x
xticklabels({'10kHz', '100kHz', '1MHz', '10MHz'}); % Define os rótulos correspondentes

% Ajustar a aparência dos rótulos do eixo x
xtickangle(45); % Rotaciona os rótulos para melhor legibilidade, se necessário

%linkaxes(h,'xy');
toc;


%save("fotorreceptor_com_filtro_azul_sem_calibracao_V2.mat","Dtt","Z","teta");
%%

[Z_corrigida,teta_corrigida] = calibracao(Z,teta,Dtt);



% Plotar os dados
% Criar uma figura
figure;

% Primeiro subplot
h(1) = subplot(2,1,1);
semilogx(Dtt, Z_corrigida, 'o-');
xlabel('Frequência (Hz)');
ylabel('Módulo Z (\Omega)');
%ylim([49 max(Z)]);
title('Módulo Z vs Frequência');
set(gca, 'XTickLabel', []); % Remover os rótulos do eixo x
grid on;

% Definir os rótulos do eixo x para mostrar os valores desejados
%xticks([10000, 100000, 1000000, 10000000]); % Define os pontos do eixo x
%xticklabels({'10kHz', '100kHz', '1MHz', '10MHz'}); % Define os rótulos correspondentes

% Ajustar a aparência dos rótulos do eixo x
xtickangle(45); % Rotaciona os rótulos para melhor legibilidade, se necessário

% Segundo subplot
h(2) = subplot(2,1,2);
semilogx(Dtt, teta_corrigida, 'o-');
xlabel('Frequência (Hz)');
%ylim([-10 max(teta)]);
ylabel('Ângulo teta (°)');
title('Ângulo teta vs Frequência');
grid on;

% Definir os rótulos do eixo x para mostrar os valores desejados
xticks([10000, 100000, 1000000, 10000000]); % Define os pontos do eixo x
xticklabels({'10kHz', '100kHz', '1MHz', '10MHz'}); % Define os rótulos correspondentes

% Ajustar a aparência dos rótulos do eixo x
xtickangle(45); % Rotaciona os rótulos para melhor legibilidade, se necessário

%linkaxes(h,'xy');

%save("fotorreceptor_com_filtro_azul__com_calibracao_V2.mat","Dtt","Z_corrigida","teta_corrigida");
%%
%Bode em Inglês

% Plotar os dados
% Criar uma figura
figure;

% Primeiro subplot
h(1) = subplot(2,1,1);
semilogx(Dtt, Z_corrigida, 'o-');
%xlabel('Frequency (Hz)');
ylabel('Module Z (\Omega)');
title('Photoreceptor - Z Bode Plot');
set(gca, 'XTickLabel', []); % Remover os rótulos do eixo x
grid on;

% Definir os rótulos do eixo x para mostrar os valores desejados
%xticks([10000, 100000, 1000000, 10000000]); % Define os pontos do eixo x
%xticklabels({'10kHz', '100kHz', '1MHz', '10MHz'}); % Define os rótulos correspondentes

% Ajustar a aparência dos rótulos do eixo x
xtickangle(45); % Rotaciona os rótulos para melhor legibilidade, se necessário

% Segundo subplot
h(2) = subplot(2,1,2);
semilogx(Dtt, teta_corrigida, 'o-');
%xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
%title('LED - Gráfico do Ângulo teta vs Frequência');
grid on;

% Definir os rótulos do eixo x para mostrar os valores desejados
xticks([10000, 100000, 1000000, 10000000]); % Define os pontos do eixo x
xticklabels({'10kHz', '100kHz', '1MHz', '10MHz'}); % Define os rótulos correspondentes

% Ajustar a aparência dos rótulos do eixo x
xtickangle(45); % Rotaciona os rótulos para melhor legibilidade, se necessário

% Ajustar a posição do segundo subplot
pos2 = get(gca, 'Position'); % Obter a posição atual
pos2(2) = pos2(2) + 0.1; % Ajustar a posição vertical (diminuir mais para mover para baixo)
set(gca, 'Position', pos2); % Aplicar a nova posição

% Ajustar a posição do rótulo do eixo X do segundo subplot
hXLabel2 = xlabel('Frequency (Hz)');
set(hXLabel2, 'Units', 'normalized');
posXLabel2 = get(hXLabel2, 'Position');
posXLabel2(2) = posXLabel2(2) + 0.15;
set(hXLabel2, 'Position', posXLabel2);
%linkaxes(h,'xy');
