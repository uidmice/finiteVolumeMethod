function myMovie = draw(r1, r2, x, E1, E2, dt, n,K1, ita1, ita2)

nImages = length(E1);
numberOfFrames = floor(nImages/n);
fig = figure;
vidHeight = 827;
vidWidth = 1120;
allTheFrames = cell(numberOfFrames,1);
allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};
allTheColorMaps = cell(numberOfFrames,1);
allTheColorMaps(:) = {zeros(256, 3)};
myMovie = struct('cdata', allTheFrames, 'colormap', allTheColorMaps);
set(gcf, 'renderer', 'zbuffer');
for i = 1:numberOfFrames
    rmax = max(max(r1((i-1)*n+1,:)),max(r2((i-1)*n+1,:)));
    subplot(2,2,1);
    cla;
    plot(x, r1((i-1)*n+1, :))
    title("|K1|>2*ita");
    ylim([0,rmax])
    
    subplot(2,2,2);
    cla;
    plot(x, r2((i-1)*n+1, :))
    title("|K1|<2*ita")
    ylim([0,rmax])
    
    subplot(2,2,3);
    cla;
    plot(linspace(0, dt*(i-1)*n, (i-1)*n+1), E1(1:(i-1)*n+1));
    s1 = sprintf("K1=%s, ita=%.1f", K1, ita1);
    title(s1);
    xlim([0,ceil(dt*i*n)])
    xlabel("t")
    ylabel("E")
    
    subplot(2,2,4);
    cla;
    plot(linspace(0, dt*(i-1)*n, (i-1)*n+1), E2(1:(i-1)*n+1));
    s2 = sprintf("K1=%s, ita=%.1f", K1, ita2);
    title(s2);
    xlim([0,ceil(dt*i*n)])
    xlabel("t")
    ylabel("E")
    
    set(fig, 'Position', [20 20 vidWidth vidHeight]);
    drawnow
    
    frame = getframe(fig);
    myMovie(i)  = frame;
end


end