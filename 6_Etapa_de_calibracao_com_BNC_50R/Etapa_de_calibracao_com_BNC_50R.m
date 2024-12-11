%PARA O BNC 50R


%Usarei esse programa para fazer a aquisição dos três valores de Z e teta
%para as aquisições V1, V2 e V3, para depois fazer a média dos três valores
%e usar essa média para remover o erro do sistema e usar para a calibração.


tic;
clc
clear



% Listar os arquivos da pasta de interesse


for V = 1:3

V_way = num2str(V); %Para puxar de forma automática o caminho do arquivo

Arquivo = "D:\dados_usuario\Documents\GEDRE\TCC_Git\Código Matlab\6_Etapa_de_calibracao_com_BNC_50R\aquisicoes_para_calibracao\50R\V" + V_way;


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


Z_mean(V,:) = Z;
teta_mean(V,:) = teta;

end

Z = (Z_mean(1,:) + Z_mean(2,:) + Z_mean(3,:))/3;  %Faço a médias das três aquisições (V1,V2 e V3)
teta = (teta_mean(1,:) + teta_mean(2,:) + teta_mean(3,:))/3; %Faço a médias das três aquisições (V1,V2 e V3)

%PARA SALVAR O VETOR, MUDAMOS OS NOMES DA VARIÁVES PARA NÃO DAR CONFLITO
%QUANDO CARREGARMOS OS VETORES PARA 50R E 25R
Pot_50R = Pot;
S_50R = S;
FP_50R = FP;
teta_50R = teta;
Z_50R = Z;
Dtt_50R = Dtt;

%save("media_R50.mat","Dtt_50R","Z_50R","teta_50R","FP_50R","S_50R","Pot_50R"); %salva os vetores salvos da média das três aquisições 
%load("media_R50.mat"); %carregar os vetores salvos da média das três aquisições

%Converto Z para a forma retangular para poder calcular o erro

% Convertendo ângulo para radianos 

radianos_50R = deg2rad(teta_50R);

% Calculando partes real e imaginária
real_50R = Z_50R .* cos(radianos_50R);
imaginaria_50R = Z_50R .* sin(radianos_50R);

Zret_50R = real_50R + (imaginaria_50R*1i);


%ERROS
erro_abs_real_50R = real_50R - 50;
erro_per_real_50R = ((real_50R - 50)/50)*100;

%Só é possível o erro absoluto da da parte imaginária, já que a referência
%é 0i
erro_abs_imag_50R = imaginaria_50R - 0;



 Amostra = 1:length(Z_50R);
% 
% 
 tab_50R = table(Amostra', Dtt_50R', Z_50R', teta_50R', erro_abs_real_50R', erro_per_real_50R', erro_abs_imag_50R', FP_50R',   'VariableNames', {'Amostras','freq', 'Z_50R', 'teta_50R', 'erro absoluto parte real', 'erro percentualparte real (%)', 'erro absoluto parte imaginária', 'FP'});
% 
% tab_mean = table({'Média'}, mean(Pot), mean(S), mean(FP), mean(teta), mean(Z), mean(Zret), mean(freq_fft), mean(freq_t), mean(Dtt), 'VariableNames', {'Amostras','Pot', 'S', 'FP', 'teta', 'Z', 'Zret', 'freq pela fft', 'freq pelo delta tempo', 'freq pelo todo tempo'});
% 
 disp(tab_50R);
% disp(tab_mean);



% Plotar os dados
% Criar uma figura
figure;



% Primeiro subplot
h(1) = subplot(2,1,1);
semilogx(Dtt_50R, Z_50R, 'o-');
xlabel('Frequência (Hz)');
ylabel('Módulo Z (\Omega)');
title('Conector BNC 50R - Gráfico do Módulo Z vs Frequência');
grid on;


% Definir os rótulos do eixo x para mostrar os valores desejados
xticks([1000, 10000, 100000, 1000000, 10000000]); % Define os pontos do eixo x
xticklabels({'1kHz', '10kHz', '100kHz', '1MHz', '10MHz'}); % Define os rótulos correspondentes

% Ajustar a aparência dos rótulos do eixo x
xtickangle(45); % Rotaciona os rótulos para melhor legibilidade, se necessário

% Segundo subplot
h(2) = subplot(2,1,2);

semilogx(Dtt_50R, teta_50R, 'o-');
xlabel('Frequência (Hz)');
ylabel('Ângulo teta (°)');
title('Conector BNC 50R - Gráfico do Ângulo teta vs Frequência');
grid on;

% Definir os rótulos do eixo x para mostrar os valores desejados
xticks([1000, 10000, 100000, 1000000, 10000000]); % Define os pontos do eixo x
xticklabels({'1kHz', '10kHz', '100kHz', '1MHz', '10MHz'}); % Define os rótulos correspondentes

% Ajustar a aparência dos rótulos do eixo x
xtickangle(45); % Rotaciona os rótulos para melhor legibilidade, se necessário

%linkaxes(h,'xy');
toc;

%%

%PARA O BNC 25R

tic;
clc
clear



% Listar os arquivos da pasta de interesse


for V = 1:3

V_way = num2str(V); %Para puxar de forma automática o caminho do arquivo

Arquivo = "D:\dados_usuario\Documents\GEDRE\TCC_Git\Código Matlab\6_Etapa_de_calibracao_com_BNC_50R\aquisicoes_para_calibracao\25R\V" + V_way;


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


Z_mean(V,:) = Z;
teta_mean(V,:) = teta;

end

Z = (Z_mean(1,:) + Z_mean(2,:) + Z_mean(3,:))/3;  %Faço a médias das três aquisições (V1,V2 e V3)
teta = (teta_mean(1,:) + teta_mean(2,:) + teta_mean(3,:))/3; %Faço a médias das três aquisições (V1,V2 e V3)

%PARA SALVAR O VETOR, MUDAMOS OS NOMES DA VARIÁVES PARA NÃO DAR CONFLITO
%QUANDO CARREGARMOS OS VETORES PARA 50R E 25R
Pot_25R = Pot;
S_25R = S;
FP_25R = FP;
teta_25R = teta;
Z_25R = Z;
Dtt_25R = Dtt;

%save("media_R25.mat","Dtt_25R","Z_25R","teta_25R","FP_25R","S_25R","Pot_25R"); %salva os vetores salvos da média das três aquisições 
load("media_R25.mat"); %carregar os vetores salvos da média das três aquisições

%Converto Z para a forma retangular para poder calcular o erro

% Convertendo ângulo para radianos 

radianos_25R = deg2rad(teta_25R);

% Calculando partes real e imaginária
real_25R = Z_25R .* cos(radianos_25R);
imaginaria_25R = Z_25R .* sin(radianos_25R);

Zret_25R = real_25R + (imaginaria_25R*1i);


%ERROS
erro_abs_real_25R = real_25R - 25;
erro_per_real_25R = ((real_25R - 25)/25)*100;

%Só é possível o erro absoluto da da parte imaginária, já que a referência
%é 0i
erro_abs_imag_25R = imaginaria_25R - 0;



 Amostra = 1:length(Z_25R);
% 
% 
 tab_25R = table(Amostra', Dtt_25R', Z_25R', teta_25R', erro_abs_real_25R', erro_per_real_25R', erro_abs_imag_25R', FP_25R',   'VariableNames', {'Amostras','freq', 'Z_25R', 'teta_25R', 'erro absoluto parte real', 'erro percentualparte real (%)', 'erro absoluto parte imaginária', 'FP'});
% 
% tab_mean = table({'Média'}, mean(Pot), mean(S), mean(FP), mean(teta), mean(Z), mean(Zret), mean(freq_fft), mean(freq_t), mean(Dtt), 'VariableNames', {'Amostras','Pot', 'S', 'FP', 'teta', 'Z', 'Zret', 'freq pela fft', 'freq pelo delta tempo', 'freq pelo todo tempo'});
% 
 disp(tab_25R);
% disp(tab_mean);



% Plotar os dados
% Criar uma figura
figure;



% Primeiro subplot
h(1) = subplot(2,1,1);
semilogx(Dtt_25R, Z_25R, 'o-');
xlabel('Frequência (Hz)');
ylabel('Módulo Z (\Omega)');
title('Conector BNC 25R - Gráfico do Módulo Z vs Frequência');
grid on;


% Definir os rótulos do eixo x para mostrar os valores desejados
xticks([1000, 10000, 100000, 1000000, 10000000]); % Define os pontos do eixo x
xticklabels({'1kHz', '10kHz', '100kHz', '1MHz', '10MHz'}); % Define os rótulos correspondentes

% Ajustar a aparência dos rótulos do eixo x
xtickangle(45); % Rotaciona os rótulos para melhor legibilidade, se necessário

% Segundo subplot
h(2) = subplot(2,1,2);

semilogx(Dtt_25R, teta_25R, 'o-');
xlabel('Frequência (Hz)');
ylabel('Ângulo teta (°)');
title('Conector BNC 25R - Gráfico do Ângulo teta vs Frequência');
grid on;

% Definir os rótulos do eixo x para mostrar os valores desejados
xticks([1000, 10000, 100000, 1000000, 10000000]); % Define os pontos do eixo x
xticklabels({'1kHz', '10kHz', '100kHz', '1MHz', '10MHz'}); % Define os rótulos correspondentes

% Ajustar a aparência dos rótulos do eixo x
xtickangle(45); % Rotaciona os rótulos para melhor legibilidade, se necessário

%linkaxes(h,'xy');
toc;

%%










