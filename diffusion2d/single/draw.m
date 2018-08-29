function draw(x,y,z)

angle={[0,0],[-90,0],[0 90],[-37.5,30]};
ttl={'front view','left view','top view', 'perspective view'};
for i=1:4
    subplot(2,2,i);
    cla;
    surf(x,y,z,'EdgeColor','none','FaceColor','interp');
    view(angle{i});title(ttl{i});
end


end