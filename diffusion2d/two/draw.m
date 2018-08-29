

function draw(x,y,z1,z2,limit)

angle={[0,0],[-90,0],[0 90],[-37.5,30]};
ttl={'front view','left view','top view', 'perspective view'};
for i=1:4
    subplot(3,4,i);
    cla;
    surf(x,y,z1,'EdgeColor','none','FaceColor','interp');
    hold on
    surf(x,y,z2,'EdgeColor','none','FaceColor','interp');
    axis(limit);
    view(angle{i});title(ttl{i});
end
for i=5:8
    subplot(3,4,i);
    surf(x,y,z1,'EdgeColor','none','FaceColor','interp');
    axis(limit);
    view(angle{i-4});
end
for i=9:12
    subplot(3,4,i);
    surf(x,y,z2,'EdgeColor','none','FaceColor','interp');
    axis(limit);
    view(angle{i-8});
end


end



