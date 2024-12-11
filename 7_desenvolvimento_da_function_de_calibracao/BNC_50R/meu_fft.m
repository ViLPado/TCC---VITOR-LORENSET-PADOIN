function [f,mod,fase] = meu_fft(sinal,Fs,ratioFMAXreduce)
%
% function [f,mod] = meu_fft(sinal,Fs,ratioFMAXreduce)
%
% Calcula a FFT e retorna a escala de frequencias e modulo 
% parametros meu_fft(sinal,Fs,ratioFMAXreduce)
% e retorna valor de amplitude(mod) dois vetores:  [f,mod]
%                                                  Hz  linear
% ratioFMAXreduce é um fator de escala que reduz a resolução da saída e
% acelera calculo
if nargin<3
    ratioFMAXreduce=1;
end

Y = fft(sinal,round(length(sinal)/ratioFMAXreduce));    
Namostras=length(Y);
%Y = fftshift(Y); % desloca a componente DC para o centro do vetor
if nargout>2
    fase = angle(Y);
    P2 = fase;
    P1 = P2(1:Namostras/2+1); % separa somente metade do espectro 
    P1(2:end-1) = P1(2:end-1); % ajusta por ter cortado o espectro
    fase = [P1];
    clear P1;
    clear P2;
end
P2 = abs(Y/Namostras);
P1 = P2(1:round(Namostras/2)+1)*sqrt(2); % separa somente metade do espectro 
P1(2:end-1) = sqrt(2)*P1(2:end-1); % ajusta por ter cortado o espectro
% para manter a energia dos componentes AC multiplica por raiz de 2
%P1(1)=0; % remove nivel DC
f = Fs*(0:round(Namostras/2))/Namostras;
mod = [P1];
clear P1;
clear P2;

end