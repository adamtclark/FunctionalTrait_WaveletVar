function plotcat(dispbin, meanmat, scale, dispnamesbin)
    colvar=size(dispbin,2);
    collst=0:(1/colvar):1;
    collst=collst(2:end);
    revcollst=collst(end:-1:1);
    for(i=1:size(dispbin,2))
        loglog(scale, quantile(meanmat(:,:,i)', [0.159, 0.841])', 'linestyle', '--', 'color', [collst(i) 0 revcollst(i)]); hold on;
        loglog(scale, quantile(meanmat(:,:,i)', [0.5])', 'marker', '.', 'linestyle', 'none', 'color', [collst(i) 0 revcollst(i)])
        tmp=meanmat(:,:,i);
        text(max(scale(:)), max(tmp(:)), dispnamesbin(i,:), 'color', [collst(i) 0 revcollst(i)], 'fontsize', 15)
    end
    rln=refline(0,1);
    set(rln, 'Color', 'k', 'LineStyle', '--', 'Marker', '.')
    xlim([2, max(scale)])
    ylim([0.8, round(max(meanmat(:)))])

    xlabel('scale, m','fontsize',15)
    ylabel('wavelet variance','fontsize',15)

    set(gca,'fontsize',15)
    set(gca,'linewidth',2)
    
    NumTicks = 5;
    L = get(gca,'XLim');
    tcks=unique(round(exp(log(L(1)):(log(L(2))-log(L(1)))/NumTicks:log(L(2)))./5).*5);
    set(gca,'XTick',tcks)
    
    L = get(gca,'YLim');
    tcks=unique(round(exp(log(L(1)):(log(L(2))-log(L(1)))/NumTicks:log(L(2)))./2).*2);
    set(gca,'YTick',tcks)
    hold off
end