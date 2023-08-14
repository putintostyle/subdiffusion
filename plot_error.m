function plot_error(alpha, t, array, label, fname, isSave)
    figure;
    for i =1:size(array,1)
        sz = 20;
        p = loglog(t(end:-1:1), array(i, :));
        p.LineStyle = "-";
        p.Marker = '.';
        p.MarkerSize = sz;
        p.DisplayName = '\alpha = '+string(alpha(i));
        hold on
%         x = linspace(t(end), t(1), 1000);
%         y = array(i, end).*(x.^alpha(i));
%         loglog(x, y, ':');
    end
    xlabel("time step")
    ylabel("error")
    grid on;
    title("Temporal error")
    legend('Location','southeast')
    if isSave
        print(gcf,append(label,'.png'),'-dpng','-r900');
        hold off;
        saveas(gcf, append(fname, '.png'));
    else
        hold off
    end
end
