function [t] = meu_tempo_osciloscopio_DHO804(arquivo1)

%     Este código serve para o osciloscópio DHO804
%     Pois ao fazer aquisições em Mesure, não há coluna de tempo.
% Apenas são apresentados o tempo inicial (t0) e o tempo entre as
% amostras (tInc)
%     Portanto, este código gera um vetor de tempo para as N amo-
% stras do arquivo CSV analisado.

%o CÓDIGO PRINCIPAL JÁ TEM O CLEAR

arquivo = arquivo1;

Dados = arquivo; %"D:\dados_usuario\Documents\GEDRE\TCC I\Código Matlab\Códigos polidos mais ou menos\2CHS 1M pontos.csv"; % Para puxar as informações da tabela 

%Verificar existência
if exist(Dados, 'file') == 2
    disp('O arquivo de CSV existe e está acessível.');
else
    disp('O arquivo de CSV não existe ou não está acessível, por favor atribua um arquivo.');
end

%%
dados = fopen(Dados, 'r');  % Abrir o arquivo para leitura
linha = textscan(dados, '%s', 1, 'Delimiter', '\n');  % Ler a primeira linha como uma célula de strings
fclose(dados);  % Fechar o arquivo

cabecalho = linha{1};  % Tem que fazer isso para extrair a primeira linha da célula, dai pode mostrar a primeira linha
                 % Só que ela sai como um string abc = {'cabeçalho'} e
                 % precisamos fazer com que seja abc = {c a b e ç a l h o}
                 % para poder encontrar o sinal de '='

abc = strjoin(cabecalho);

%%

virgula = find(abc==',');
igual = find(abc=='=');


if strcmp(abc(6),'C')==0 %Para saber se estou trabalhando com dois canais ou apenas um canal do osciloscópio
    disp('CH1V');
    t0 = abc(igual(1)+1:virgula(2)-1);
    t1 = abc(igual(2)+1:virgula(3)-1);

    t0 = str2num(t0);
    t1 = str2num(t1); %se não aparecer é porque é muito pequeno
else
    disp('CH1V e CH2V');
    t0 = abc(igual(1)+1:virgula(3)-1);
    t1 = abc(igual(2)+1:virgula(4)-1);

    t0 = str2num(t0);
    t1 = str2num(t1);
end
%%

VIdados1 = readtable(Dados);

t = zeros(1,height(VIdados1));
t(1)=t0;

for i=2:1:height(VIdados1)
    t(i)=t1+t(i-1);
end
