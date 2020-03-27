function plot_fn(N,t,c,cp,plot_ct)

mark=["d","s","o"];
line=["-",":","-."];
cole=["r","b","g"];

    figure
    for i=1:N
        hold on
        plot(t, c(:,i),mark(i)+cole(i),...
             t,cp(:,i),line(i)+cole(i))
    end


end