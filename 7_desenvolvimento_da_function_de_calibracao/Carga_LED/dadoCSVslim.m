function [Pot,S,FP,teta,Z,Zret,freq_fundamental,f_fundamental,Dtt,Num_ciclos] = dadoCSVslim(local_arquivo,n)

% A objetivo deste código é:
% - Importar os dados de V e I de um arquivo CSV
% - A partir do Vrms e do Irms é cálculado o módulo da impedância
%   e a paritr do FP = P/S é possível
%   obtermos o ângulo da impedância (cos^-1 FP)
% OBS: isso para apenas uma frequência
%
% O arquivo CSV foi baseado no modelo que é fornecido pelo oscilo-
% cópio DHO804 (caso seja de outra forma, o códido deve ser modificado).
% O arquivo CSV deve conter:
% - Primeira linha: CH1V,CH2V,t0 =XXXXe-XX, tInc = XXXXe-XX
% - CH1V: Tensão (V)
% - CH2V: Corrente (A)
%
% O formato da função é:
% [Pot,S,FP,teta,Z,Zret,freq] = dadoCSVslim(local_arquivo,n)
%
% ENTRADAS:
% local_arquivo -> Caminho para o arquivo CSV
% n -> Números de linhas inicias que devem ser ignoradas do arquivo CSV
% Para o DHO804 n=1
%
% SAÍDAS:
% Pot -> Potência Ativa (W)
% S -> Potência Aparente (VA)
% FP -> Fator de Potência
% teta -> Ângulo teta (Ver se está em graus)
% Z -> Impedância na forma polar (Ohm)
% Zret -> Impedância na forma retangular (Ohm)
% freq -> Frequência dos sinais analisados (ver se está em Hz)

VI = local_arquivo; %'D:\dados_usuario\Documents\GEDRE\TCC I\Código Matlab\Códigos polidos mais ou menos\Amostras\Com carga cap mais res\Carga 100k0.csv'; %edit(V) -> para ver as informações da tabela 

% Usando a função exist para verificar a existência do arquivo
if exist(VI, 'file') == 2
    disp('O arquivo de tensão e corrente existe e está acessível.');
else
    disp('O arquivo de tensão e corrente não existe ou não está acessível, por favor atribua um arquivo.');
end

%%

%
% VER O SINAL DO OSCILOSCÓPIO
%

VIdados_all = readtable(VI, "NumHeaderLines",n); %"NumHeaderLines",16 -> Para dizer a paritr de qual linha da tabela CSV do osciloscóio que começa os números (informações) de verdade
                                                 %open Vdados_all  -> Para ver a matriz de dados ou edit(V) -> para ver as informações da tabela
VIdados1 = table2array(VIdados_all); %Passa de uma tabela para um array, mais fácil de manipular

%keyboard


[t] = meu_tempo_osciloscopio_DHO804(VI);
Vamp = zeros(1,height(VIdados_all)); %Criar vetor como o número de espeços necessários agiliza o processo
Iamp = zeros(1,height(VIdados_all));

i = 0;

for i = 1:1:height(VIdados1)   %Passar de uma coluna para um array   
     Vamp(1,i)= VIdados1(i,1);
     Iamp(1,i)= VIdados1(i,2);
end

% 
% plot(t*10^3, Vamp,'b', t*10^3, Iamp);
% title('Dados do Osciloscópio de Tensão e Corrente');
% xlabel('Tempo (ms)');
% ylabel('Amplitude'); 
%%
%Potência Ativa: multiplica VxI ponto a ponto, soma e divide pelo número de
%amostras (média)

%Potência Aparente: Vrms x Irms, para achar o rms, é a variável ao quadrado
%(Ex: V -> V^2), e faz raiz da média, Vrms = sqrt(mean(V.^2)). OU só Vrms =
%rms(V)

%Negócio de plotar depois que passar por zero (lembrar que as amostras tem
%que ter ciclos inteiros)

%%
%
% VER A PASSAGEM POR ZERO
%


% figure
% hold on
% plot(diff(Vamp>0)>0,'r') %para plotar um sinal em cima do outro, posso clicar só em um gráfico e depois plotar o outro)

%%
%
% FILTRAGEM
%

% Frequência fundamental do sinal
taxa_amostragem = 1/((t(end)-t(1))/(height(VIdados_all)-1)); %Intervalo entre as amostras (s)
[freq,mod]=meu_fft(Vamp,taxa_amostragem);

%mod(1) = 0; %removendo intêncinalmente o nível DC
position_mod_max = find(mod==max(mod));
freq_fundamental = freq(position_mod_max);

% Parâmetros do filtro
frequencia_corte = freq_fundamental*3;  % Frequência de corte em Hertz
ordem = 1000;  % Ordem do filtro                                                 %Aumentar a ordem do filtro gera defasagem entre o sinal original e o sinal filtrado


% Gerando um sinal de exemplo
taxa_amostragem = 1/((t(end)-t(1))/(height(VIdados_all)-1));  % Taxa de amostragem em Hz, quanto tempo leva cada amostra, por exemplo, o osciloscópio de 2GHz, pode ter um tempo entre as amostras de até 1/2GHz
tempo = t;
V_sinal_original = Vamp;
I_sinal_original = Iamp;

% Criando o filtro passa-baixas
filtro_passa_baixas = fir1(ordem, frequencia_corte/(taxa_amostragem/2), 'low');

% Aplicando o filtro ao sinal
V_sinal_filtrado = filter(filtro_passa_baixas, 1, V_sinal_original);
I_sinal_filtrado = filter(filtro_passa_baixas, 1, I_sinal_original);

% Exibindo o sinal original e o sinal filtrado
% figure;
% h(1) = subplot(2,1,1);
% plot(tempo, V_sinal_original, tempo, I_sinal_original);
% title('Sinal Original');
% xlabel('Tempo (s)');
% ylabel('Amplitude');
% 
% h(2) = subplot(2,1,2);
% plot(tempo, V_sinal_filtrado, 'g', tempo, I_sinal_filtrado);
% title('Sinal Filtrado (Passa-Baixas)');
% xlabel('Tempo (s)');
% ylabel('Amplitude');
% 
% linkaxes(h,'xy')

%%
% figure
% plot(V_sinal_filtrado)
% hold on
% plot(I_sinal_filtrado)
% hold on
% plot(diff(V_sinal_filtrado>0)>0,'r') 

%%

x=diff(V_sinal_filtrado>0)>0;
%plot(diff(sinal_filtrado>0)>0,'r');

% plot(V_sinal_filtrado)
% hold on
% plot(Vamp)
% hold on
% plot(I_sinal_filtrado)
% hold on
% plot(Iamp)
% hold on
% plot(x,'r')
%%
%plot(cumsum(x)) %quantas vezes passa por zero, depois da parea implementar um jeito de usar isso para ver se está passando apenas uma vez por zero por ciclo,
                %achar um jeito de calcular a frequência da fundamental
                %(fft talvez), e saber ou definir um número X de ciclos, e
                %aprtir disso, saber quantas vezes tem que passar por zero.
                %Se for passar mais de uma vez, da para automatizar para
                %aumentar a ordem do filtro.
                %OU da para calcular a diferença de tempo entre cada ciclo
                %e achar a frequência, se quiser.
                %

%%

%
% AJUSTES PARA RETIRAR OS EFEITOS DO FILTRO NOS PONTOS INICAIS, que ocorrem de 0:ordem
%

t1 = t(ordem:end);
V_sinal_filtrado1 = V_sinal_filtrado(ordem:end); %Tirar o efeito do filtro nas primeiras amostras
I_sinal_filtrado1 = I_sinal_filtrado(ordem:end); %Tirar o efeito do filtro nas primeiras amostras
Vamp1 = Vamp(ordem:end);
Iamp1 = Iamp(ordem:end);
%%

x=diff(V_sinal_filtrado>0)>0;
x1 = x(ordem-1:end);

% plot(t1, V_sinal_filtrado1)
% hold on
% %plot(t1,Vamp1)
% plot(t1, I_sinal_filtrado1)
% hold on
% %plot(t1,Iamp1)
% hold on
% plot(t1,x1,'r')

%%

%
% ENCONTRAR A PRIMEIRA PASSAGEM POR ZERO
%

pz = find(x1); %passagens por zero (vindo do negativo para o positivo)

%Primeira passagem por zero

pz(1);

%
% FREQUÊNCIA FUNDAMNETAL
%

t1(pz(1));
t1(pz(2));

Dt = t1(pz(2))-t1(pz(1));
f_fundamental = 1/Dt;

Dtt = 1/((t1(pz(end))-t1(pz(1)))/(length(pz)-1));

Num_ciclos=(length(pz)-1);
% f_fundamental_arredondada = ceil(1/(Dt*1e3));  %Em Hz, rredonda para cima, mas arredoda na faixa dos kHz, Ex: se 4490 Hz -> 5000 Hz

%%

%
% NÚMERO DE CICLOS QUE QUERO USAR
%

t2 = t1(pz(1):pz(end));
V_sinal_filtrado2 = V_sinal_filtrado1(pz(1):pz(end));
Vamp2 = Vamp1(pz(1):pz(end));
I_sinal_filtrado2 = I_sinal_filtrado1(pz(1):pz(end));
Iamp2 = Iamp1(pz(1):pz(end));
x2 = x1(pz(1):pz(end));

% plot(t2, V_sinal_filtrado2)
% % hold on
% % plot(t2,Vamp2)
% hold on
% plot(t2, I_sinal_filtrado2)
% % hold on
% % plot(t2,Iamp2)
% hold on
% plot(t2,x2,'r')

%%

%
% CÁLCULO DAS POTÊNCIAS
%

% ATIVA
P = (V_sinal_filtrado2.*I_sinal_filtrado2); %multiplicação ponto a ponto dos sinais de V e I filtrados e condicionados para um número fechado de ciclos
Pot = mean(P);

%disp(['Potência Ativa: ' num2str(Pot) ' W']);

% APARENTE
Vrms = rms(V_sinal_filtrado2);
Irms = rms(I_sinal_filtrado2);

S = Vrms*Irms;

%disp(['Potência Aparente: ' num2str(S) ' VA']);

%%

%
% CÁLCULA DA IMPEDÂNCIA
%



FP = Pot/S;
teta = acos(round(FP, 10))*180/pi;

if (V_sinal_filtrado2(1)-I_sinal_filtrado2(1)) < 0 %-0.000000000001   %I adiantada -> Cap
    teta = -teta;
else
    teta = teta;
end

%disp(['Ângulo teta: ' num2str(teta, '%.4f') ' °']);


Z = Vrms/Irms;

%disp(['Z_Polar: (' num2str(Z, '%.4f') ',' num2str(teta, '%.4f') '°)']);

% Convertendo ângulo para radianos
radianos = deg2rad(teta);

% Calculando partes real e imaginária
real = Z * cos(radianos);
imaginaria = Z * sin(radianos);

Zret = real + (imaginaria*1i);

%disp(['Z_Retangular: ' num2str(real) ' + j*' num2str(imaginaria) char(937)]);

%disp(['Frequência fundamental ftt: ' num2str(freq_fundamental/1000, '%.4f') 'kHz']);
%disp(['Frequência fundamental tempo: ' num2str(f_fundamental/1000, '%.4f') 'kHz']);
%disp(['Frequência fundamental: ' num2str(Dtt/1000, '%.8f') 'kHz']);

%%
% plot(meu_fft(Vamp,taxa_amostragem))
% semilogy(freq,mod)
% plot(freq,mod)
% [freq,mod]=meu_fft(Vamp,taxa_amostragem); % Ver um jeito de pegar o maior módulo e a paritr disso descobrir a freq, depois disso, usar umas 3*freq como a freq de corte
%
% help addpath % Colocar a do meu_fft