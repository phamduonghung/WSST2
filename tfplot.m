function h=tfplot(Ts,typ,freq,time)

N = size(Ts,2);
if nargin<4
    time = (0:N-1)/N;
end
if nargin<3
    if typ=='lin' 
        freq = 0:N/2-1;
    else
        freq = 1:N/2;
    end
end
        

switch typ
    case 'lin'
        % Linear TF plot
        h=figure();
        imagesc((0:N-1)/N,0:N/2-1,log(1+abs(Ts)));
        set(gca,'YDir','normal');
%         colormap('normal');
%         set(gcf,'colormap',flipud(get(gcf,'colormap')));
        set(gca, 'xtick', []) ;
        xlabel('t');ylabel('\eta');
    case 'log'
        h = figure();ha = axes();
        % Log frequencies
        flin = 1:length(freq);
        imagesc(time,flin,log(0.01+(abs(Ts))),'Parent',ha);
        set(ha,'YDir','normal');
        %set(gca, 'xtick', []) ;
        %xlabel('t');
        ylabel('1/a (log)'); 
        %colormap(1-gray);
        %set(gcf,'colormap',flipud(get(gcf,'colormap')));
        % log ticks
        set(ha,'YTick',(linspace(flin(1),flin(end),6)),'YTickLabel',round(freq(floor(linspace(1,length(freq),6))),2));
end